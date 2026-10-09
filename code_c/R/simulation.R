# Table 11 extension. True parameters are used only in artificial-data experiments.
c2_simulation_signature <- function(cfg,source_identity=c2_source_identity()) {
  c2_md5(list(version="simulation-v1",cfg=cfg,source=source_identity))
}

c2_simulation_scenarios <- function() {
  list(list(name="mfBm",kappa=NULL),list(name="mfOU_slow",kappa=c(.2,.4)),
       list(name="mfOU_fast",kappa=c(2,4)))
}

c2_simulation_draw <- function(R,N,mu,seed,brownian=FALSE) {
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit(if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv))
  set.seed(seed)
  Y <- matrix(drop(t(R)%*%rnorm(nrow(R))),ncol=length(mu),byrow=TRUE)
  if(brownian) Y <- rbind(rep(0,length(mu)),Y)
  if(nrow(Y)!=N) stop("simulation covariance/draw dimensions disagree")
  Y <- sweep(Y,2,mu,"+");colnames(Y) <- paste0("asset",seq_along(mu));Y
}

c2_simulation_parameter_rows <- function(p,dgp,replicate,model,truth,mu=NULL) {
  rows <- list()
  add <- function(asset,parameter,value,true_value) {
    rows[[length(rows)+1L]] <<- data.frame(dgp=dgp,replicate=replicate,model=model,
      asset=asset,parameter=parameter,estimate=value,true_value=true_value,
      error=value-true_value,truth_scope=if(parameter %in% c("kappa","mu"))
        if(is.na(true_value)) "not_defined_for_mfBm" else "stationary_mfOU_parameter" else
        "driving_mfBm_parameter",stringsAsFactors=FALSE)
  }
  for(i in seq_along(p$H)) {
    add(paste0("asset",i),"H",p$H[i],truth$H[i])
    add(paste0("asset",i),"sigma",p$sigma[i],truth$sigma[i])
    if(!is.null(p$kappa)) add(paste0("asset",i),"kappa",p$kappa[i],
      if(is.null(truth$kappa)) NA_real_ else truth$kappa[i])
    if(!is.null(mu)) add(paste0("asset",i),"mu",mu[i],
      if(is.null(truth$kappa)) NA_real_ else truth$mu[i])
  }
  add("asset1:asset2","rho",p$rho[1,2],truth$rho[1,2])
  c2_bind(rows)
}

c2_run_simulation <- function(cfg,cache_dir=NULL,progress=FALSE) {
  n <- as.integer(cfg$sim_n); B <- as.integer(cfg$sim_reps); hs <- as.integer(cfg$horizons)
  delta <- cfg$sim_delta
  if(length(n)!=1L || is.na(n) || n<65L || length(B)!=1L || is.na(B) || B<1L ||
     !length(hs) || any(is.na(hs)|hs<1L) || anyDuplicated(hs) ||
     length(delta)!=1L || !is.finite(delta) || delta<=0 || !is.finite(cfg$seed))
    stop("invalid simulation configuration: need >=65 training levels, positive replications and horizons")
  q <- if(is.null(cfg$quadrature)) 128L else cfg$quadrature
  qfit <- if(is.null(cfg$fit_quadrature)) 64L else cfg$fit_quadrature
  signature <- c2_simulation_signature(cfg)
  if(!is.null(cache_dir)) {
    cache_dir <- file.path(cache_dir,signature)
    dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  }
  N <- n+max(hs); times <- (0:(N-1L))*delta
  truth <- list(H=c(.2,.4),sigma=c(1,1),rho=matrix(c(1,.4,.4,1),2),mu=c(-2,-2))
  bp <- mfbm_params(truth$H,truth$rho,0,truth$sigma)
  bw <- mfbm_forecast_weights(bp,(1:(n-1L))*delta,(n-1L+hs)*delta,target=1)
  scenarios <- c2_simulation_scenarios();cells <- list();hits <- 0L
  for(scenario in scenarios) {
    stationary <- !is.null(scenario$kappa)
    tp <- truth;tp$kappa <- scenario$kappa
    op <- if(stationary) c2_mfou_params(tp$H,tp$kappa,tp$sigma,tp$rho,tp$mu) else NULL
    # Shared exact level covariance and Cholesky across replications of this scenario.
    S <- if(stationary) c2_mfou_cov_matrix(op,times,q) else mfbm_level_cov_matrix(bp,times[-1L])
    R <- chol(S)
    for(rep in seq_len(B)) {
      path <- if(is.null(cache_dir)) NULL else file.path(cache_dir,sprintf("%s_%06d.rds",scenario$name,rep))
      if(!is.null(path) && file.exists(path)) {
        cell <- readRDS(path)
        if(!identical(cell$signature,signature) || !identical(cell$dgp,scenario$name) || !identical(cell$replicate,rep))
          stop("invalid simulation checkpoint identity")
        hits <- hits+1L
      } else {
        seed <- strtoi(substr(c2_md5(list(seed=cfg$seed,scenario=scenario$name,replicate=rep)),1,7),16L)
        Y <- c2_simulation_draw(R,N,tp$mu,seed,brownian=!stationary)
        train <- Y[seq_len(n),,drop=FALSE]
        actual_log <- Y[n+hs,1L];actual <- exp(actual_log)
        param_rows <- list();diagnostic_rows <- list();failure_rows <- list();forecast_rows <- list()
        models <- c("mfBm_known_driver","mfBm_estimated","mfOU_known","mfOU_estimated")
        for(model in models) {
          known <- model %in% c("mfBm_known_driver","mfOU_known")
          not_applicable <- model=="mfOU_known" && !stationary
          result <- NULL;reason <- "";status <- "ok"
          fitted_par <- NULL;fit_diagnostics <- NULL
          if(not_applicable) {
            status <- "not_applicable";reason <- "mfBm DGP has no positive true kappa for a stationary mfOU oracle"
          } else {
            result <- tryCatch({
              if(model=="mfBm_known_driver") {
                Z <- sweep(train[-1L,,drop=FALSE],2,train[1L,])
                m <- train[1L,1L]+drop(crossprod(as.vector(t(Z)),bw$W))
                list(forecast=exp(m+bw$msfe/2),log_mean=m,log_variance=bw$msfe)
              } else if(model=="mfBm_estimated") {
                p <- window_mfbm_params_calendar(train,delta)
                fitted_par <- p
                w <- mfbm_forecast_weights(p,(1:(n-1L))*delta,(n-1L+hs)*delta,target=1)
                Z <- sweep(train[-1L,,drop=FALSE],2,train[1L,])
                m <- train[1L,1L]+drop(crossprod(as.vector(t(Z)),w$W))
                # The published B forecaster remains the authoritative RV baseline.
                fc <- forecast_mfbm_window(train,hs,delta)
                list(forecast=fc,log_mean=m,log_variance=w$msfe,par=p)
              } else if(model=="mfOU_known") {
                c2_mfou_condition(train,hs,delta,op,estimate_mean=FALSE,quadrature=q)
              } else {
                fit <- c2_fit_mfou(train,delta,quadrature=qfit)
                fitted_par <- fit$par;fit_diagnostics <- fit$diagnostics
                out <- c2_mfou_condition(train,hs,delta,fit$par,estimate_mean=TRUE,quadrature=q)
                out$par <- fit$par;out$diagnostics <- fit$diagnostics;out
              }
            },error=function(e) {reason <<- conditionMessage(e);NULL})
            if(is.null(result)) status <- "failed"
            if(!is.null(result) && (any(!is.finite(c(result$forecast,result$log_mean,result$log_variance))) ||
                any(result$forecast<=0))) {
              status <- "failed";reason <- "nonfinite or nonpositive simulation forecast";result <- NULL
            }
          }
          if(!known && !is.null(fitted_par))
            param_rows[[length(param_rows)+1L]] <- c2_simulation_parameter_rows(fitted_par,scenario$name,rep,model,tp,
              mu=if(model=="mfOU_estimated") result$mu else NULL)
          if(!is.null(fit_diagnostics)) {
            d <- fit_diagnostics;d$dgp <- scenario$name;d$replicate <- rep;d$model <- model
            diagnostic_rows[[length(diagnostic_rows)+1L]] <- d
          }
          if(status=="failed") failure_rows[[length(failure_rows)+1L]] <- data.frame(
            dgp=scenario$name,replicate=rep,model=model,reason=reason,stringsAsFactors=FALSE)
          f <- if(is.null(result)) rep(NA_real_,length(hs)) else result$forecast
          m <- if(is.null(result)) rep(NA_real_,length(hs)) else result$log_mean
          v <- if(is.null(result)) rep(NA_real_,length(hs)) else result$log_variance
          specification <- if(model=="mfBm_known_driver" && stationary) "misspecified_covariance_known_driver" else
            if(grepl("^mfBm",model) && stationary) "misspecified_mfBm" else
              if(grepl("^mfOU",model) && !stationary) "misspecified_stationary_mfOU" else "correct_model_family"
          forecast_rows[[length(forecast_rows)+1L]] <- data.frame(dgp=scenario$name,replicate=rep,seed=seed,
            model=model,parameter_knowledge=if(known) "known_simulation_only" else "estimated_training_only",
            model_specification=specification,h=hs,forecast=f,forecast_log=m,
            model_log_variance=v,actual=actual,actual_log=actual_log,error=f-actual,log_error=m-actual_log,
            status=status,reason=reason,stringsAsFactors=FALSE)
        }
        cell <- list(signature=signature,dgp=scenario$name,replicate=rep,
          forecasts=c2_bind(forecast_rows),parameter_errors=c2_bind(param_rows),
          diagnostics=c2_bind(diagnostic_rows),failures=c2_bind(failure_rows))
        if(!is.null(path)) {
          partial <- tempfile("cell_",tmpdir=dirname(path),fileext=".partial")
          saveRDS(cell,partial,version=3)
          if(!file.rename(partial,path)) stop("atomic simulation checkpoint rename failed")
        }
      }
      cells[[length(cells)+1L]] <- cell
      if(progress) cat(sprintf("simulation %s replication %d/%d%s\n",scenario$name,rep,B,
        if(!is.null(path) && hits>0L) " (checkpoint enabled)" else ""))
    }
  }
  pull <- function(name) c2_bind(lapply(cells,`[[`,name))
  forecasts <- pull("forecasts");parameters <- pull("parameter_errors")
  failures <- pull("failures")
  if(!ncol(failures)) failures <- data.frame(dgp=character(),replicate=integer(),
    model=character(),reason=character(),stringsAsFactors=FALSE)
  summaries <- c2_simulation_metrics(forecasts,parameters,requested_reps=B)
  c(list(forecasts=forecasts,parameter_errors=parameters,diagnostics=pull("diagnostics"),
    failures=failures,cfg=cfg,signature=signature,n_cache_hits=hits),summaries)
}

c2_simulation_metrics <- function(forecasts,parameter_errors=NULL,requested_reps=NULL) {
  f <- forecasts
  needed <- c("dgp","replicate","h","model","status","actual","forecast","log_error","error")
  if(!all(needed %in% names(f))) stop("missing simulation forecast columns")
  key <- paste(f$dgp,f$replicate,f$h,f$model,sep="::")
  if(anyDuplicated(key)) stop("duplicate simulation forecast keys")
  mcse <- function(x) if(length(x)>1L) sd(x)/sqrt(length(x)) else NA_real_
  losses <- function(x) {
    ratio <- x$actual/x$forecast
    list(MSFE=x$error^2,MSFE_log=x$log_error^2,QLIKE=ratio-log(ratio)-1)
  }
  requested <- function(x) if(is.null(requested_reps)) length(unique(x$replicate)) else as.integer(requested_reps)
  groups <- split(f,paste(f$dgp,f$model,f$h,sep="::"));rows <- list()
  for(g in groups) {
    ok <- g$status=="ok" & is.finite(g$forecast) & g$forecast>0 & is.finite(g$actual) & g$actual>0 &
      is.finite(g$error) & is.finite(g$log_error)
    L <- losses(g[ok,,drop=FALSE]);ns <- sum(ok)
    avg <- function(x) if(length(x)) mean(x) else NA_real_
    rows[[length(rows)+1L]] <- data.frame(dgp=g$dgp[1],model=g$model[1],h=g$h[1],
      n_requested=requested(g),n_success=ns,n_failed=sum(g$status=="failed"),
      n_not_applicable=sum(g$status=="not_applicable"),
      MSFE=avg(L$MSFE),RMSFE=sqrt(avg(L$MSFE)),MCSE_MSFE=mcse(L$MSFE),
      MSFE_log=avg(L$MSFE_log),RMSFE_log=sqrt(avg(L$MSFE_log)),MCSE_MSFE_log=mcse(L$MSFE_log),
      QLIKE=avg(L$QLIKE),MCSE_QLIKE=mcse(L$QLIKE),stringsAsFactors=FALSE)
  }
  pairs <- list(
    mfBm_estimation_error=c("mfBm_known_driver","mfBm_estimated"),
    mfOU_estimation_error=c("mfOU_known","mfOU_estimated"),
    estimated_model_comparison=c("mfBm_estimated","mfOU_estimated"),
    known_model_comparison=c("mfBm_known_driver","mfOU_known"))
  comparisons <- list()
  for(dgp in unique(f$dgp)) for(h in unique(f$h)) for(comparison in names(pairs)) {
    p <- pairs[[comparison]];g <- f[f$dgp==dgp & f$h==h,,drop=FALSE]
    ref <- g[g$model==p[1] & g$status=="ok",,drop=FALSE]
    cand <- g[g$model==p[2] & g$status=="ok",,drop=FALSE]
    reps <- intersect(ref$replicate,cand$replicate)
    a <- ref[match(reps,ref$replicate),,drop=FALSE];b <- cand[match(reps,cand$replicate),,drop=FALSE]
    if(length(reps) && !isTRUE(all.equal(a$actual,b$actual,tolerance=0))) stop("simulation paired outcomes disagree")
    La <- losses(a);Lb <- losses(b)
    for(loss in names(La)) {
      delta_loss <- Lb[[loss]]-La[[loss]]
      base <- if(length(reps)) mean(La[[loss]]) else NA_real_
      gain <- if(length(reps) && is.finite(base) && base>0) 100*(1-mean(Lb[[loss]])/base) else NA_real_
      comparisons[[length(comparisons)+1L]] <- data.frame(dgp=dgp,h=h,comparison=comparison,
        reference=p[1],candidate=p[2],loss=loss,n_requested=requested(g),n_paired=length(reps),
        mean_loss_difference=if(length(reps)) mean(delta_loss) else NA_real_,
        MCSE_loss_difference=mcse(delta_loss),improvement_pct=gain,stringsAsFactors=FALSE)
    }
  }
  parameter_summary <- list()
  if(!is.null(parameter_errors) && nrow(parameter_errors)) {
    ps <- split(parameter_errors,paste(parameter_errors$dgp,parameter_errors$model,
      parameter_errors$asset,parameter_errors$parameter,sep="::"))
    for(p in ps) {
      good <- is.finite(p$error);e <- p$error[good]
      parameter_summary[[length(parameter_summary)+1L]] <- data.frame(dgp=p$dgp[1],model=p$model[1],
        asset=p$asset[1],parameter=p$parameter[1],truth_scope=p$truth_scope[1],
        n_requested=if(is.null(requested_reps)) length(unique(p$replicate)) else requested_reps,
        n_estimated=nrow(p),n_with_truth=length(e),true_value=p$true_value[1],
        bias=if(length(e)) mean(e) else NA_real_,RMSE=if(length(e)) sqrt(mean(e^2)) else NA_real_,
        MCSE_bias=mcse(e),stringsAsFactors=FALSE)
    }
  }
  list(metrics=c2_bind(rows),comparisons=c2_bind(comparisons),parameter_summary=c2_bind(parameter_summary))
}
