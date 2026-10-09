# Original C question retained; predictors and eligibility come from reviewed B.
c2_selection_scores <- function(Y,target,candidates,delta,min_coverage) {
  rows <- list()
  for(asset in candidates) {
    pair <- Y[,c(target,asset),drop=FALSE]
    eligible <- code_b_eligible_origins(pair,nrow(pair),nrow(pair),min_coverage)
    if(!eligible$eligible) next
    fit <- tryCatch(window_mfbm_params_calendar(pair,delta),error=function(e) NULL)
    if(is.null(fit)) next
    rho <- fit$rho[1,2]; gap <- abs(diff(fit$H))
    rows[[length(rows)+1L]] <- data.frame(asset=asset,rho=rho,H_target=fit$H[1],H_asset=fit$H[2],
      H_gap=gap,correlation_score=abs(rho),theory_score=abs(rho)*gap)
  }
  c2_bind(rows)
}

c2_selection_cell <- function(panel,cfg,target,t) {
  w <- cfg$window; k <- cfg$selection_k
  Y <- as.matrix(panel[(t-w+1L):t,cfg$selection_universe,drop=FALSE]);Y[!is.finite(Y)] <- NA_real_
  meta <- data.frame(target=target,origin_date=panel$date[t],origin_index=t,window=w)
  skip <- function(reason) list(forecasts=data.frame(),selected_assets=data.frame(),audit=cbind(meta,reason=reason))
  scores <- c2_selection_scores(Y,target,setdiff(colnames(Y),target),cfg$delta,cfg$min_window_coverage)
  fixed <- cfg$selection_fixed[[target]]
  if(nrow(scores)<k || !all(fixed %in% scores$asset)) return(skip("fixed_or_candidate_ineligible"))
  sets <- list(fixed=fixed,
    correlation=scores$asset[order(-scores$correlation_score,scores$asset)][seq_len(k)],
    theory_guided=scores$asset[order(-scores$theory_score,scores$asset)][seq_len(k)])
  # Pairwise candidate eligibility alone does not imply eligibility of a joint set.
  eligible <- vapply(sets,function(s) code_b_eligible_origins(Y[,c(target,s),drop=FALSE],w,w,cfg$min_window_coverage)$eligible,TRUE)
  if(!all(eligible)) return(skip("selected_joint_set_ineligible"))
  keep <- t+cfg$horizons<=nrow(panel)
  ix <- which(keep)
  keep[ix] <- is.finite(panel[[target]][t+cfg$horizons[ix]]) &
    panel$date[t+cfg$horizons[ix]]>=cfg$start & panel$date[t+cfg$horizons[ix]]<=cfg$end
  if(!any(keep)) return(skip("no_evaluable_target"))
  predictions <- lapply(sets,function(s) forecast_mfbm_calendar(Y[,c(target,s),drop=FALSE],cfg$horizons,cfg$delta))
  predictions$fBm <- forecast_mfbm_calendar(Y[,target,drop=FALSE],cfg$horizons,cfg$delta)
  past <- as.matrix(panel[1:t,target,drop=FALSE])
  predictions$HAR <- forecast_har_window(past,cfg$horizons,w,"har",na_avg=TRUE)
  fc <- c2_bind(lapply(names(predictions),function(rule) cbind(meta,strategy=rule,
    h=cfg$horizons[keep],target_date=panel$date[t+cfg$horizons[keep]],calendar="trading",
    forecast=predictions[[rule]][keep],actual=exp(panel[[target]][t+cfg$horizons[keep]]))))
  if(any(!is.finite(fc$forecast)|fc$forecast<=0)) stop("invalid selection forecast")
  selected <- c2_bind(lapply(names(sets),function(rule) cbind(meta,strategy=rule,rank=seq_len(k),scores[match(sets[[rule]],scores$asset),])))
  list(forecasts=fc,selected_assets=selected,audit=cbind(meta,reason="paired_realized",n_candidate_eligible=nrow(scores)))
}

c2_run_selection <- function(panel,cfg,origins=NULL,cache_dir=NULL,progress=FALSE) {
  c2_validate_panel(panel,cfg$selection_universe)
  if(is.null(origins)) origins <- c2_origins(panel,cfg)
  signature <- c2_md5(list(panel=panel[,c("date",cfg$selection_universe)],cfg=cfg,source=c2_source_identity()))
  if(!is.null(cache_dir)) {cache_dir<-file.path(cache_dir,signature);dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)}
  rows <- list();hits<-0L
  for(tk in cfg$selection_targets) for(t in origins) {
    path <- if(!is.null(cache_dir)) file.path(cache_dir,paste0(tk,"_",t,".rds")) else NULL
    if(!is.null(path)&&file.exists(path)) {cell<-readRDS(path);hits<-hits+1L} else {
      cell<-c2_selection_cell(panel,cfg,tk,t)
      if(!is.null(path)) {tmp<-paste0(path,".partial");saveRDS(cell,tmp);if(!file.rename(tmp,path)) stop("checkpoint rename failed")}
    }
    rows[[length(rows)+1L]]<-cell
    if(progress) cat("selection",tk,as.character(panel$date[t]),nrow(cell$forecasts),"rows\n")
  }
  # Some audit rows have a candidate count, some skips do not.
  audits<-lapply(rows,function(z) {if(!"n_candidate_eligible" %in% names(z$audit)) z$audit$n_candidate_eligible<-NA_integer_;z$audit})
  list(forecasts=c2_bind(lapply(rows,`[[`,"forecasts")),selected_assets=c2_bind(lapply(rows,`[[`,"selected_assets")),
    audit=c2_bind(audits),cfg=cfg,signature=signature,n_cache_hits=hits)
}

c2_selection_metrics <- function(fc) {
  if(!nrow(fc)) return(data.frame())
  groups <- split(fc,paste(fc$target,fc$h))
  c2_bind(lapply(groups,function(g) {
    base <- g[g$strategy=="fixed",]
    c2_bind(lapply(unique(g$strategy),function(rule) {
      x <- g[g$strategy==rule,];idx<-match(c2_forecast_key(x),c2_forecast_key(base))
      if(anyNA(idx) || any(x$actual!=base$actual[idx])) stop("selection comparison keys differ")
      msfe<-mean(loss_se(x$actual,x$forecast));qlike<-mean(loss_qlike(x$actual,x$forecast))
      data.frame(target=x$target[1],h=x$h[1],strategy=rule,n=nrow(x),msfe=msfe,rmsfe=sqrt(msfe),qlike=qlike,
        msfe_gain_pct=100*(1-msfe/mean(loss_se(base$actual[idx],base$forecast[idx]))),
        qlike_gain_pct=100*(1-qlike/mean(loss_qlike(base$actual[idx],base$forecast[idx]))),scope="descriptive_selection_heuristic")
    }))
  }))
}
