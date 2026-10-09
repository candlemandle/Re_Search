# Prespecified uncertainty sensitivity on an already executed run; never refits models.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<1L) stop("Usage: Rscript code_c/scripts/compare_blocks.R results/code_c_v2/dense_study/mfou/empirical_run.rds")
options(code_c.root=normalizePath("code_c"))
for(f in sort(list.files("code_c/R",pattern="[.]R$",full.names=TRUE))) source(f)
ab<-c2_load_ab(if(length(args)>1L) args[2] else NULL)
x<-readRDS(args[1]);out<-dirname(args[1]);rows<-list()
for(block in c(20L,40L,60L)) {
  z<-c2_compare_models(x$forecasts,x$cfg$mode,x$cfg$bootstrap_reps,block,x$cfg$min_bootstrap_blocks,x$cfg$seed)
  z$block_requested<-block;rows[[length(rows)+1L]]<-z
}
c2_write_csv(c2_bind(rows),file.path(out,"bootstrap_block_sensitivity.csv"))
cat("Saved paired pointwise block sensitivity; unavailable cells remain unavailable.\n")
