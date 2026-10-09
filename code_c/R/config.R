# Additive entry point: A/B remain owned by colleagues.
c2_ab_root <- function(root=NULL) {
  if(is.null(root)) root <- Sys.getenv("CODE_C_AB_ROOT","")
  if(!nzchar(root)) root <- if(dir.exists("part B version 2/R")) "part B version 2" else "."
  normalizePath(root,mustWork=TRUE)
}

c2_load_ab <- function(root=NULL) {
  root <- c2_ab_root(root)
  source(file.path(root,"R/code_b_config.R")); source_code_b(root)
  required <- c("forecast_mfbm_calendar","window_mfbm_params_calendar","code_b_eligible_origins","code_b_calendar")
  if(!all(vapply(required,exists,TRUE,mode="function")) ||
      !identical(code_b_config("full")$calendar,"trading"))
    stop("Code C v2 needs reviewed B v2 with its trading-calendar API")
  root
}

c2_config <- function(mode=c("smoke","study","full")) {
  mode <- match.arg(mode); b <- code_b_config("full")
  assets <- b$dj30_assets
  list(mode=mode,seed=20261009L,delta=b$delta,window=b$window,
    horizons=b$horizons,assets=assets,sets=lapply(seq_along(assets),function(k) assets[1:k]),
    sensitivity_windows=c(250L,1000L),calendar="trading",min_window_coverage=b$min_window_coverage,
    start=as.Date(b$dj30_periods$full[1]),end=as.Date(b$dj30_periods$full[2]),
    smoke_origins=3L,study_start=as.Date("2019-01-02"),study_end=as.Date("2021-12-31"),study_step=21L,
    sim_n=if(mode=="smoke") 80L else 500L,
    sim_reps=switch(mode,smoke=4L,study=30L,full=200L),sim_delta=1/250,
    bootstrap_reps=if(mode=="full") 5000L else 2000L,bootstrap_block_days=20L,
    min_bootstrap_blocks=10L,quadrature=128L,fit_quadrature=64L,
    selection_targets=c("AAPL","MSFT","JPM"),selection_universe=names(code_a_config("full")$dj30),
    selection_k=4L,selection_fixed=list(AAPL=c("ALD","AMGN","AXP","BA"),
      MSFT=c("AAPL","ALD","AMGN","AXP"),JPM=c("AAPL","ALD","AMGN","AXP")))
}

c2_bind <- function(rows) {
  if(!length(rows)) return(data.frame())
  x <- do.call(rbind,rows); rownames(x) <- NULL; x
}

c2_md5 <- function(x) {
  f <- tempfile(fileext=".rds"); on.exit(unlink(f)); saveRDS(x,f,version=3)
  unname(tools::md5sum(f))
}

c2_source_identity <- function(ab_root=getOption("code_b.source_root",".")) {
  ab <- sort(list.files(file.path(ab_root,"R"),pattern="[.]R$",full.names=TRUE))
  own <- sort(list.files(file.path(getOption("code_c.root","code_c"),"R"),pattern="[.]R$",full.names=TRUE))
  if(!length(ab)||!length(own)) stop("cannot fingerprint A/B and C sources")
  setNames(unname(tools::md5sum(c(ab,own))),c(paste0("AB/",basename(ab)),paste0("C/",basename(own))))
}

c2_origins <- function(panel,cfg) {
  o <- seq.int(cfg$window,nrow(panel)-1L)
  o <- o[panel$date[o]>=cfg$start-max(cfg$horizons)*2 & panel$date[o]<cfg$end]
  if(cfg$mode=="smoke") o <- o[unique(round(seq(1,length(o),length.out=cfg$smoke_origins)))]
  if(cfg$mode=="study") {
    o <- o[panel$date[o]>=cfg$study_start & panel$date[o]<=cfg$study_end]
    o <- o[seq.int(1L,length(o),by=cfg$study_step)]
  }
  o
}

c2_scenarios <- function(cfg) {
  out <- lapply(cfg$sets,function(a) list(assets=a,window=cfg$window))
  for(w in cfg$sensitivity_windows) for(k in intersect(c(1L,5L),lengths(cfg$sets)))
    out[[length(out)+1L]] <- list(assets=cfg$assets[1:k],window=w)
  out
}
