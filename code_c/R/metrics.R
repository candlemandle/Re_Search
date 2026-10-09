c2_forecast_key <- function(x) paste(x$origin_date,x$target_date,x$h)

c2_pairs <- function(x) {
  a <- x[x$model=="mfBm",]; b <- x[x$model=="mfOU",]
  ka <- c2_forecast_key(a); kb <- c2_forecast_key(b)
  if(anyDuplicated(ka)||anyDuplicated(kb)) stop("duplicate forecast key")
  ia <- which(ka %in% kb); ib <- match(ka[ia],kb)
  a <- a[ia,,drop=FALSE]; b <- b[ib,,drop=FALSE]
  if(any(a$actual!=b$actual)) stop("paired actuals differ")
  list(base=a,proposed=b)
}

c2_metrics <- function(fc) {
  if(!nrow(fc)) return(data.frame())
  fc$period <- "full"
  extra <- list(fc)
  for(p in c("period1","period2")) {
    inside <- if(p=="period1") fc$target_date<=as.Date("2017-04-11") else
      fc$target_date>=as.Date("2017-04-12") & fc$target_date<=as.Date("2021-07-30")
    z <- fc[inside,]; if(!nrow(z)) next
    z$period <- p;extra[[length(extra)+1L]] <- z
  }
  z <- c2_bind(extra)
  groups <- split(z,paste(z$dimension,z$window,z$h,z$period))
  c2_bind(lapply(groups,function(g) {
    pair <- c2_pairs(g)
    c2_bind(lapply(c("mfBm","mfOU"),function(m) {
      x <- if(m=="mfBm") pair$base else pair$proposed
      if(!nrow(x)) return(data.frame())
      data.frame(dimension=x$dimension[1],window=x$window[1],h=x$h[1],period=x$period[1],model=m,n=nrow(x),
        msfe=mean(loss_se(x$actual,x$forecast)),rmsfe=sqrt(mean(loss_se(x$actual,x$forecast))),
        qlike=mean(loss_qlike(x$actual,x$forecast)),first_target=min(x$target_date),last_target=max(x$target_date))
    }))
  }))
}

c2_window_comparison <- function(fc,reference_window=500L) {
  if(!nrow(fc)) return(data.frame())
  out <- list()
  for(d in unique(fc$dimension)) for(m in unique(fc$model)) for(h in unique(fc$h)) {
    x <- fc[fc$dimension==d & fc$model==m & fc$h==h,]
    a <- x[x$window==reference_window,]
    for(w in setdiff(unique(x$window),reference_window)) {
      b <- x[x$window==w,];ka<-c2_forecast_key(a);kb<-c2_forecast_key(b)
      ia<-which(ka %in% kb);ib<-match(ka[ia],kb)
      if(!length(ia)) next
      if(any(a$actual[ia]!=b$actual[ib])) stop("window actuals differ")
      for(metric in c("MSFE","QLIKE")) {
        f<-if(metric=="MSFE") loss_se else loss_qlike
        bl<-mean(f(a$actual[ia],a$forecast[ia]));pl<-mean(f(b$actual[ib],b$forecast[ib]))
        out[[length(out)+1L]]<-data.frame(dimension=d,model=m,h=h,reference_window=reference_window,
          alternative_window=w,metric=metric,n=length(ia),reference_loss=bl,alternative_loss=pl,
          gain_pct=100*(1-pl/bl),scope="matched_keys_descriptive_window_sensitivity")
      }
    }
  }
  c2_bind(out)
}

c2_loss_paths <- function(fc) {
  if(!nrow(fc)) return(data.frame())
  c2_bind(lapply(split(fc,paste(fc$dimension,fc$window,fc$h)),function(x) {
    z<-c2_pairs(x);a<-z$base;b<-z$proposed;if(!nrow(a)) return(data.frame())
    o<-order(a$origin_date);a<-a[o,];b<-b[o,]
    data.frame(dimension=a$dimension,window=a$window,h=a$h,origin_date=a$origin_date,target_date=a$target_date,
      cumulative_se_difference=cumsum(loss_se(a$actual,a$forecast)-loss_se(b$actual,b$forecast)),
      cumulative_qlike_difference=cumsum(loss_qlike(a$actual,a$forecast)-loss_qlike(b$actual,b$forecast)))
  }))
}

c2_compare_models <- function(fc,mode="smoke",B=2000L,block_days=20L,min_blocks=10L,seed=20261009L,period="full") {
  if(!nrow(fc)) return(data.frame())
  groups <- split(fc,paste(fc$dimension,fc$window,fc$h))
  c2_bind(lapply(groups,function(g) {
    pair <- c2_pairs(g); a <- pair$base; b <- pair$proposed
    if(!nrow(a)) return(data.frame())
    order <- order(a$origin_date); a <- a[order,];b <- b[order,]; n <- nrow(a)
    gaps <- if("origin_index" %in% names(a)) diff(a$origin_index) else as.numeric(diff(a$origin_date))
    step <- if(length(gaps)) median(gaps) else NA_real_
    has_trading_index<-"origin_index" %in% names(a)
    positions<-if(has_trading_index) a$origin_index-min(a$origin_index)+1L else seq_len(n)
    n_calendar<-max(positions);coverage<-n/n_calendar
    # Keep missing pairs on the calendar instead of compressing nonconsecutive days.
    # At least 90% paired coverage is an a priori guard; sparse study grids cannot infer.
    dense <- has_trading_index && length(gaps)>0 && step==1 && coverage>=.9
    block <- max(block_days,a$h[1]); available <- mode!="smoke" && dense && n>=block*min_blocks
    intervals <- matrix(NA_real_,2,2)
    if(available) {
      set.seed(seed+a$dimension[1]*1000L+a$window[1]+a$h[1])
      base_losses <- proposed_losses <- matrix(NA_real_,n_calendar,2)
      base_losses[positions,]<-cbind(loss_se(a$actual,a$forecast),loss_qlike(a$actual,a$forecast))
      proposed_losses[positions,]<-cbind(loss_se(b$actual,b$forecast),loss_qlike(b$actual,b$forecast))
      draws <- replicate(B,{idx <- block_bootstrap_index(n_calendar,block)
        100*(1-colMeans(proposed_losses[idx,,drop=FALSE],na.rm=TRUE)/colMeans(base_losses[idx,,drop=FALSE],na.rm=TRUE))})
      good<-colSums(is.finite(draws))==2L
      if(sum(good)<.95*B) stop("too few finite paired calendar bootstrap draws")
      intervals <- apply(draws[,good,drop=FALSE],1,quantile,probs=c(.025,.975))
    }
    c2_bind(lapply(seq_len(2),function(j) {
      f <- if(j==1L) loss_se else loss_qlike
      bl <- mean(f(a$actual,a$forecast)); pl <- mean(f(b$actual,b$forecast))
      data.frame(dimension=a$dimension[1],window=a$window[1],h=a$h[1],period=period,metric=c("MSFE","QLIKE")[j],n=n,
        baseline_loss=bl,proposed_loss=pl,gain_pct=100*(1-pl/bl),ci_low=intervals[1,j],ci_high=intervals[2,j],
        uncertainty_status=if(mode=="smoke") "unavailable_smoke" else if(available) "exploratory_pointwise" else "unavailable_sparse_or_insufficient_blocks",
        bootstrap_B=if(available) B else 0L,block_trading_days=block,origin_step_median=step,
        n_calendar_origins=n_calendar,paired_calendar_coverage=coverage,
        bootstrap_calendar=if(has_trading_index) "trading_positions_missing_pairs_retained" else "no_trading_index_inference_unavailable",
        n_origin_gaps=if(length(gaps)) sum(gaps!=1) else 0L,inferential_scope=if(available) "paired_exploratory_no_multiple_test_control" else "descriptive_only")
    }))
  }))
}

c2_compare_scopes <- function(fc,cfg) {
  out<-list()
  for(period in c("full","period1","period2")) {
    keep<-switch(period,full=rep(TRUE,nrow(fc)),period1=fc$target_date<=as.Date("2017-04-11"),
      period2=fc$target_date>=as.Date("2017-04-12") & fc$target_date<=as.Date("2021-07-30"))
    x<-fc[keep,,drop=FALSE];if(!nrow(x)) next
    out[[length(out)+1L]]<-c2_compare_models(x,cfg$mode,cfg$bootstrap_reps,cfg$bootstrap_block_days,
      cfg$min_bootstrap_blocks,cfg$seed,period)
  }
  c2_bind(out)
}
