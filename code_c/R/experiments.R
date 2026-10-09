c2_validate_panel <- function(panel,assets) {
  if(!inherits(panel$date,"Date") || anyNA(panel$date) || anyDuplicated(panel$date) ||
      is.unsorted(panel$date,strictly=TRUE) || !all(assets %in% names(panel))) stop("invalid input panel")
  invisible(TRUE)
}

c2_comparison_assets <- function(cfg,window) {
  sets <- Filter(function(sc) sc$window==window,c2_scenarios(cfg))
  unique(unlist(lapply(sets,`[[`,"assets"),use.names=FALSE))
}

c2_joint_eligibility <- function(panel,cfg,window,origins) {
  L <- as.matrix(panel[,c2_comparison_assets(cfg,window),drop=FALSE])
  L[!is.finite(L)] <- NA_real_
  code_b_eligible_origins(L,origins,window,cfg$min_window_coverage)
}

c2_empirical_cell <- function(panel,cfg,scenario,t) {
  a <- scenario$assets; w <- scenario$window; d <- length(a)
  meta <- data.frame(dimension=d,window=w,origin_date=panel$date[t],origin_index=t)
  audit <- function(reason,h=NA_integer_) cbind(meta,h=h,reason=reason)
  L <- as.matrix(panel[,a,drop=FALSE]); L[!is.finite(L)] <- NA_real_
  eligible <- c2_joint_eligibility(panel,cfg,w,t)
  if(!eligible$eligible) return(list(audit=audit(eligible$reason),forecasts=data.frame(),parameters=data.frame(),failures=data.frame()))
  keep <- t+cfg$horizons<=nrow(panel)
  reason <- rep("outside_panel",length(keep))
  reason[keep] <- ifelse(is.finite(L[t+cfg$horizons[keep],1]),"realized","missing_actual")
  valid <- which(keep)
  keep[valid] <- is.finite(L[t+cfg$horizons[valid],1]) & panel$date[t+cfg$horizons[valid]]>=cfg$start & panel$date[t+cfg$horizons[valid]]<=cfg$end
  reason[valid[!keep[valid] & reason[valid]!="missing_actual"]] <- "outside_evaluation_period"
  audits <- list(); if(any(!keep)) audits[[1]] <- audit(reason[!keep],cfg$horizons[!keep])
  if(!any(keep)) return(list(audit=c2_bind(audits),forecasts=data.frame(),parameters=data.frame(),failures=data.frame()))
  Y <- L[(t-w+1L):t,,drop=FALSE]
  baseline <- forecast_mfbm_calendar(Y,cfg$horizons,cfg$delta)
  t0 <- proc.time()[[3]]
  proposed <- tryCatch({
    fit <- c2_fit_mfou(Y,cfg$delta,quadrature=cfg$fit_quadrature)
    pred <- c2_mfou_condition(Y,cfg$horizons,cfg$delta,fit$par,quadrature=cfg$quadrature)
    list(fit=fit,pred=pred)
  },error=function(e) list(error=conditionMessage(e)))
  elapsed <- proc.time()[[3]]-t0
  row <- function(name,pred,log_mean=NA_real_,log_variance=NA_real_) {
    if(any(!is.finite(pred)|pred<=0)) stop("nonfinite/nonpositive ",name," forecast")
    cbind(meta,model=name,calendar="trading",h=cfg$horizons[keep],
      target_date=panel$date[t+cfg$horizons[keep]],train_start=panel$date[t-w+1L],
      n_train_days=sum(complete.cases(Y)),forecast=pred[keep],actual=exp(L[t+cfg$horizons[keep],1]),
      log_mean=rep(log_mean,length.out=length(keep))[keep],log_variance=rep(log_variance,length.out=length(keep))[keep])
  }
  fc <- list(row("mfBm",baseline)); parameters <- failures <- data.frame()
  if(!is.null(proposed$error)) {
    failures <- cbind(meta,model="mfOU",reason=proposed$error,elapsed_seconds=elapsed)
    audits[[length(audits)+1L]] <- audit("mfOU_estimation_or_numerical_failure")
  } else {
    p <- proposed$pred
    fc[[2]] <- row("mfOU",p$forecast,p$log_mean,p$log_variance)
    parameters <- cbind(meta,proposed$fit$diagnostics,mu=p$mu,
      cholesky_diagonal_ratio=p$cholesky_diagonal_ratio,elapsed_seconds=elapsed)
    audits[[length(audits)+1L]] <- audit("paired_realized",cfg$horizons[keep])
  }
  list(forecasts=c2_bind(fc),parameters=parameters,failures=failures,audit=c2_bind(audits))
}

c2_run_empirical <- function(panel,cfg,origins=NULL,cache_dir=NULL,progress=FALSE) {
  c2_validate_panel(panel,cfg$assets)
  if(is.null(origins)) origins <- c2_origins(panel,cfg)
  if(!length(origins) || anyNA(origins) || anyDuplicated(origins) || any(origins<cfg$window|origins>=nrow(panel)))
    stop("invalid empirical origins")
  signature <- c2_md5(list(panel=panel[,c("date",cfg$assets)],cfg=cfg,source=c2_source_identity()))
  if(!is.null(cache_dir)) {cache_dir <- file.path(cache_dir,signature);dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)}
  rows <- list(); hits <- 0L
  for(sc in c2_scenarios(cfg)) for(t in origins) {
    id <- paste0("w",sc$window,"_d",length(sc$assets),"_o",t)
    path <- if(!is.null(cache_dir)) file.path(cache_dir,paste0(id,".rds")) else NULL
    if(!is.null(path) && file.exists(path)) {cell <- readRDS(path);hits <- hits+1L} else {
      cell <- c2_empirical_cell(panel,cfg,sc,t)
      if(!is.null(path)) {
        tmp <- paste0(path,".partial");saveRDS(cell,tmp);if(!file.rename(tmp,path)) stop("checkpoint rename failed")
      }
    }
    rows[[length(rows)+1L]] <- cell
    if(progress) cat(id,if(!is.null(path)&&file.exists(path)) "saved/available" else "computed",";",nrow(cell$forecasts),"rows\n")
  }
  pull <- function(name) c2_bind(lapply(rows,`[[`,name))
  list(forecasts=pull("forecasts"),parameters=pull("parameters"),audit=pull("audit"),
    failures=pull("failures"),cfg=cfg,origins=origins,signature=signature,n_cache_hits=hits,
    protocol_version="joint_trading_v2")
}

# Past-only calendar alignment of previously computed cells. Retained forecasts
# are unchanged; the reanalysis script preserves the original RDS before use.
c2_apply_joint_calendar <- function(x,panel) {
  if(identical(x$protocol_version,"joint_trading_v2")) return(x)
  removed <- 0L
  for(w in unique(vapply(c2_scenarios(x$cfg),`[[`,integer(1),"window"))) {
    e <- c2_joint_eligibility(panel,x$cfg,w,x$origins)
    excluded <- x$origins[!e$eligible]
    if(!length(excluded)) next
    for(name in c("forecasts","parameters","failures","audit")) {
      z <- x[[name]]
      if(!nrow(z)) next
      keep <- !(z$window==w & z$origin_index %in% excluded)
      if(name=="forecasts") removed <- removed+sum(!keep)
      x[[name]] <- z[keep,,drop=FALSE]
    }
    for(sc in Filter(function(s) s$window==w,c2_scenarios(x$cfg))) {
      z <- data.frame(dimension=length(sc$assets),window=w,
        origin_date=panel$date[excluded],origin_index=excluded,h=NA_integer_,reason=e$reason[!e$eligible])
      x$audit <- c2_bind(list(x$audit,z))
    }
  }
  x$protocol_version <- "joint_trading_v2"
  x$calendar_alignment <- list(removed_forecast_rows=removed,rule="past-only shared eligibility per training window")
  x
}
