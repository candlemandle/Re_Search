# From repository root: Rscript code_c/scripts/run_experiments.R smoke [all|mfou|selection|simulation]
args<-commandArgs(trailingOnly=TRUE)
mode<-if(length(args)) match.arg(args[1],c("smoke","study","full")) else "smoke"
scope<-if(length(args)>=2 && !startsWith(args[2],"--")) match.arg(args[2],c("all","mfou","selection","simulation")) else "all"
option<-function(name,default=NULL) {v<-args[startsWith(args,paste0("--",name,"="))];if(length(v)>1) stop("duplicate option ",name);if(length(v)) substring(v,nchar(name)+4) else default}
# One BLAS worker gives reproducible practical timings across dense Gaussian solves.
if(!nzchar(Sys.getenv("OPENBLAS_NUM_THREADS"))) {
  status<-system2(file.path(R.home("bin"),"Rscript"),shQuote(c("code_c/scripts/run_experiments.R",args)),
    env=c("OPENBLAS_NUM_THREADS=1","OMP_NUM_THREADS=1","VECLIB_MAXIMUM_THREADS=1"))
  quit(save="no",status=status)
}
options(code_c.root=normalizePath("code_c"))
for(f in sort(list.files("code_c/R",pattern="[.]R$",full.names=TRUE))) source(f)
ab_root<-c2_load_ab(option("ab-root"))
cfg<-c2_config(mode)
if(!is.null(option("dimensions"))) {
  dims<-as.integer(strsplit(option("dimensions"),",",fixed=TRUE)[[1]])
  if(!length(dims)||anyNA(dims)||anyDuplicated(dims)||any(!dims %in% seq_along(cfg$assets))) stop("invalid dimensions")
  cfg$sets<-lapply(sort(dims),function(d) cfg$assets[seq_len(d)])
}
if(!is.null(option("sensitivity"))) {
  value<-option("sensitivity")
  windows<-if(value=="none") integer() else as.integer(strsplit(value,",",fixed=TRUE)[[1]])
  if(anyNA(windows)||anyDuplicated(windows)||any(!windows %in% c(250L,1000L))) stop("invalid sensitivity windows")
  cfg$sensitivity_windows<-windows
}
if(!is.null(option("study-step"))) cfg$study_step<-as.integer(option("study-step"))
if(!is.null(option("study-start"))) cfg$study_start<-as.Date(option("study-start"))
if(!is.null(option("study-end"))) cfg$study_end<-as.Date(option("study-end"))
if(is.na(cfg$study_step)||cfg$study_step<1L||anyNA(c(cfg$study_start,cfg$study_end))||cfg$study_start>cfg$study_end)
  stop("invalid prespecified study calendar")
panel<-load_logvol_panel("dj30",file.path(ab_root,"data/processed"))
out<-option("output",file.path("results/code_c_v2",mode));dir.create(out,recursive=TRUE,showWarnings=FALSE)
writeLines(c(paste("AB root:",ab_root),paste("Mode:",mode),paste("Scope:",scope),capture.output(sessionInfo())),file.path(out,"session_info.txt"))
id<-c2_source_identity(ab_root);c2_write_csv(data.frame(path=names(id),md5=unname(id)),file.path(out,"source_identity.csv"))
c2_write_csv(data.frame(path="data/processed/dj30_logvol.csv",md5=unname(tools::md5sum(file.path(ab_root,"data/processed/dj30_logvol.csv")))),file.path(out,"data_identity.csv"))
saveRDS(cfg,file.path(out,"configuration.rds"))
plan<-c2_workload_plan(panel,cfg);c2_write_csv(plan,file.path(out,"workload_plan.csv"));print(plan)
if("--plan" %in% args) quit(save="no",status=0)
origins<-c2_origins(panel,cfg)
limit<-option("max-origins")
if(!is.null(limit)) {limit<-as.integer(limit);if(is.na(limit)||limit<1L) stop("invalid max-origins");origins<-head(origins,limit)}
stages<-list()
stage<-function(name,fun) {
  t0<-proc.time()[[3]];fun()
  stages[[length(stages)+1L]]<<-data.frame(stage=name,status="completed",elapsed_seconds=proc.time()[[3]]-t0)
  c2_write_csv(c2_bind(stages),file.path(out,"run_manifest.csv"))
}
if(scope %in% c("all","mfou")) stage("mfou",function() {
  x<-c2_run_empirical(panel,cfg,origins,file.path(out,"checkpoints/empirical"),progress=TRUE)
  c2_write_empirical(x,file.path(out,"mfou"))
  if(!nrow(x$parameters)) stop("No mfOU fit succeeded; inspect failure outputs")
})
if(scope %in% c("all","selection")) stage("selection",function() {
  x<-c2_run_selection(panel,cfg,origins,file.path(out,"checkpoints/selection"),progress=TRUE)
  d<-file.path(out,"selection");dir.create(d,recursive=TRUE,showWarnings=FALSE);saveRDS(x,file.path(d,"selection_run.rds"))
  for(n in c("forecasts","selected_assets","audit")) c2_write_csv(x[[n]],file.path(d,paste0(n,".csv")))
  c2_write_csv(c2_selection_metrics(x$forecasts),file.path(d,"metrics.csv"))
})
if(scope %in% c("all","simulation")) stage("simulation",function() {
  x<-c2_run_simulation(cfg,file.path(out,"checkpoints/simulation"),progress=TRUE)
  d<-file.path(out,"simulation");dir.create(d,recursive=TRUE,showWarnings=FALSE);saveRDS(x,file.path(d,"simulation_run.rds"))
  for(n in c("forecasts","parameter_errors","diagnostics","failures","metrics","comparisons","parameter_summary"))
    c2_write_csv(x[[n]],file.path(d,paste0(n,".csv")))
})
cat("Completed",mode,scope,"in",out,"\n")
