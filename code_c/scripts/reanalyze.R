# Analysis-only refresh of saved raw forecasts. Preserves original forecast provenance.
args<-commandArgs(trailingOnly=TRUE)
if(!length(args)) stop("Supply an executed empirical_run.rds path")
options(code_c.root=normalizePath("code_c"))
for(f in sort(list.files("code_c/R",pattern="[.]R$",full.names=TRUE))) source(f)
ab<-c2_load_ab(if(length(args)>1L) args[2] else NULL)
x<-readRDS(args[1])
if(!identical(x$protocol_version,"joint_trading_v2")) {
  backup<-paste0(args[1],".original_per_set.rds")
  if(!file.exists(backup) && !file.copy(args[1],backup)) stop("Cannot preserve original raw forecasts")
  panel<-load_logvol_panel("dj30",file.path(ab,"data/processed"))
  x<-c2_apply_joint_calendar(x,panel)
  cat("Aligned past-only joint calendar; removed",x$calendar_alignment$removed_forecast_rows,
      "forecast rows. Original RDS preserved at",backup,"\n")
}
c2_write_empirical(x,dirname(args[1]))
cat("Saved analysis and separate analysis_source_identity.csv; forecasts were not refitted.\n")
