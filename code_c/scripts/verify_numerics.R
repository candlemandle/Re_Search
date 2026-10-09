# Compare numerical precision and authoritative B forecasts on real 500-day windows.
args<-commandArgs(trailingOnly=TRUE)
if(!length(args)) stop("Supply the executed empirical_run.rds path")
options(code_c.root=normalizePath("code_c"))
for(f in sort(list.files("code_c/R",pattern="[.]R$",full.names=TRUE))) source(f)
ab<-c2_load_ab(if(length(args)>1L) args[2] else NULL)
x<-readRDS(args[1]);cfg<-x$cfg
panel<-load_logvol_panel("dj30",file.path(ab,"data/processed"));rows<-list()
for(d in intersect(c(1L,2L,5L),unique(x$forecasts$dimension))) {
  base<-x$forecasts[x$forecasts$dimension==d & x$forecasts$window==500L & x$forecasts$model=="mfBm",]
  if(!nrow(base)) next
  origin<-base$origin_date[1];t<-match(origin,panel$date);Y<-as.matrix(panel[(t-499L):t,cfg$assets[seq_len(d)],drop=FALSE])
  expected<-forecast_mfbm_calendar(Y,cfg$horizons,cfg$delta)
  saved<-base[base$origin_date==origin,];diffB<-max(abs(saved$forecast-expected[match(saved$h,cfg$horizons)]))
  fit<-c2_fit_mfou(Y,cfg$delta,quadrature=cfg$fit_quadrature)
  a<-c2_mfou_condition(Y,cfg$horizons,cfg$delta,fit$par,quadrature=128L)
  b<-c2_mfou_condition(Y,cfg$horizons,cfg$delta,fit$par,quadrature=256L)
  scale<-max(abs(a$forecast));relative<-max(abs(a$forecast-b$forecast))/scale
  savedOU<-x$forecasts[x$forecasts$dimension==d & x$forecasts$window==500L &
    x$forecasts$model=="mfOU" & x$forecasts$origin_date==origin,]
  if(!nrow(savedOU)) stop("No saved mfOU forecast for precision audit")
  refit_difference<-max(abs(savedOU$forecast-a$forecast[match(savedOU$h,cfg$horizons)]))/scale
  # A priori numerical tolerances, independent from observed forecast losses.
  if(diffB>1e-12 || relative>1e-4 || refit_difference>1e-8) stop("Numerical verification failed for dimension ",d)
  rows[[length(rows)+1L]]<-data.frame(dimension=d,window=500L,origin_date=origin,
    baseline_B_max_abs_difference=diffB,mfOU_128_vs_256_relative_difference=relative,
    saved_mfOU_vs_current_refit_relative_difference=refit_difference,
    baseline_tolerance=1e-12,quadrature_tolerance=1e-4,refit_tolerance=1e-8,status="passed")
}
c2_write_csv(c2_bind(rows),file.path(dirname(args[1]),"numerical_precision.csv"))
print(c2_bind(rows))
