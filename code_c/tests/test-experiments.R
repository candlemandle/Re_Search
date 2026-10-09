c2_test_panel <- function(n=140L) {
  set.seed(42)
  Y <- apply(matrix(rnorm(n*5,sd=.07),n,5),2,cumsum)-2
  data.frame(date=as.Date("2010-01-01")+seq_len(n),Y,check.names=FALSE)
}

test_that("mfBm baseline is the reviewed B trading predictor, even with a gap", {
  panel <- c2_test_panel(); names(panel)[2:6] <- c("AAPL","ALD","AMGN","AXP","BA")
  panel$ALD[40] <- NA_real_
  cfg <- c2_config("smoke"); cfg$window <- 80L; cfg$horizons <- c(1L,5L)
  cfg$sets <- list(c("AAPL","ALD")); cfg$sensitivity_windows <- integer()
  x <- c2_run_empirical(panel,cfg,origins=100L)
  expected <- forecast_mfbm_calendar(as.matrix(panel[21:100,c("AAPL","ALD")]),cfg$horizons,cfg$delta)
  fc <- x$forecasts[x$forecasts$model=="mfBm",]
  expect_equal(fc$forecast,expected,tolerance=1e-12)
  expect_equal(fc$target_date,panel$date[100L+cfg$horizons])
  expect_true(all(fc$calendar=="trading"))
  future <- panel; future$ALD[101:140] <- NA_real_; future$AAPL[110:140] <- 9
  y <- c2_run_empirical(future,cfg,origins=100L)
  expect_equal(x$forecasts,y$forecasts,tolerance=1e-12)
  numeric_parameters <- setdiff(names(x$parameters),"elapsed_seconds")
  expect_equal(x$parameters[,numeric_parameters],y$parameters[,numeric_parameters],tolerance=1e-12)
})

test_that("missing targets are dropped jointly, failures stay visible", {
  panel <- c2_test_panel(); names(panel)[2:6] <- c("AAPL","ALD","AMGN","AXP","BA")
  cfg <- c2_config("smoke"); cfg$window <- 80L; cfg$horizons <- c(1L,5L)
  cfg$sets <- list("AAPL"); cfg$sensitivity_windows <- integer()
  panel$AAPL[105] <- NA_real_
  x <- c2_run_empirical(panel,cfg,origins=100L)
  expect_true(all(x$forecasts$h==1L))
  expect_true(any(x$audit$reason=="missing_actual"))
  expect_true(all(x$forecasts$n_train_days==80L))
})

test_that("checkpoint resume is exact and changed data invalidates identity", {
  panel <- c2_test_panel(); names(panel)[2:6] <- c("AAPL","ALD","AMGN","AXP","BA")
  cfg <- c2_config("smoke"); cfg$window <- 80L; cfg$horizons <- 1L
  cfg$sets <- list("AAPL"); cfg$sensitivity_windows <- integer()
  d <- tempfile(); dir.create(d); on.exit(unlink(d,recursive=TRUE))
  x <- c2_run_empirical(panel,cfg,100L,cache_dir=d)
  y <- c2_run_empirical(panel,cfg,100L,cache_dir=d)
  expect_equal(x$forecasts,y$forecasts,tolerance=0)
  expect_equal(y$n_cache_hits,1L)
  changed <- panel; changed$AAPL[99] <- changed$AAPL[99]+.01
  z <- c2_run_empirical(changed,cfg,100L,cache_dir=d)
  expect_equal(z$n_cache_hits,0L)
})

test_that("smoke uncertainty is unavailable and losses pair exact keys", {
  fc <- data.frame(model=rep(c("mfBm","mfOU"),each=3),dimension=1L,window=500L,
    h=1L,origin_date=rep(as.Date("2020-01-01")+0:2,2),
    target_date=rep(as.Date("2020-01-01")+1:3,2),actual=.2,
    forecast=c(.19,.22,.21,.2,.21,.19))
  x <- c2_compare_models(fc,mode="smoke")
  expect_true(all(is.na(x$ci_low)&is.na(x$ci_high)))
  expect_true(all(x$uncertainty_status=="unavailable_smoke"))
  expect_equal(x$n,c(3L,3L))
  expect_equal(x$gain_pct[1],100*(1-mean((.2-c(.2,.21,.19))^2)/mean((.2-c(.19,.22,.21))^2)))
})

test_that("reports handle a missing paper subperiod without fabricating rows", {
  fc<-data.frame(model=c("mfBm","mfOU"),dimension=1L,window=500L,h=1L,
    origin_date=as.Date("2010-01-01"),target_date=as.Date("2010-01-02"),actual=.2,forecast=c(.19,.21))
  x<-c2_metrics(fc)
  expect_setequal(unique(x$period),c("full","period1"))
  expect_true(all(x$n==1L))
})

test_that("CLI plan honors prespecified dimensions and paths with spaces", {
  root<-dirname(getOption("code_c.root"));previous<-setwd(root);on.exit(setwd(previous))
  out<-file.path(tempdir(),"C plan with spaces");on.exit(unlink(out,recursive=TRUE),add=TRUE)
  args<-c("code_c/scripts/run_experiments.R","study","mfou","--plan","--dimensions=1,2", "--sensitivity=none",
    paste0("--ab-root=",getOption("code_b.source_root")),paste0("--output=",out))
  log<-tempfile();on.exit(unlink(log),add=TRUE)
  status<-system2(file.path(R.home("bin"),"Rscript"),shQuote(args),stdout=log,stderr=log)
  expect_equal(status,0L,info=paste(readLines(log),collapse="\n"))
  plan<-read.csv(file.path(out,"workload_plan.csv"))
  expect_setequal(plan$dimension,c(1L,2L));expect_true(all(plan$window==500L))
})

test_that("daily bootstrap preserves missing pair positions on the trading grid", {
  day<-rep((1:80)[-31],2);n<-79L
  fc<-data.frame(model=rep(c("mfBm","mfOU"),each=n),dimension=2L,window=500L,h=5L,
    origin_index=day,origin_date=as.Date("2020-01-01")+day,target_date=as.Date("2020-01-01")+day+5,
    actual=.2,forecast=c(.2+.01*sin(1:n),.2+.012*cos(1:n)))
  x<-c2_compare_models(fc,"study",B=100L,block_days=5L,min_blocks=5L)
  expect_true(all(x$uncertainty_status=="exploratory_pointwise"))
  expect_true(all(is.finite(x$ci_low)&is.finite(x$ci_high)))
  expect_true(all(x$n_calendar_origins==80L))
  expect_equal(x$paired_calendar_coverage,rep(79/80,2),tolerance=1e-12)
  expect_true(all(x$n_origin_gaps==1L))
})

test_that("all failed mfOU fits still produce a usable failure report", {
  fc<-data.frame(model="mfBm",dimension=1L,window=500L,h=1L,
    origin_date=as.Date("2010-01-01"),target_date=as.Date("2010-01-02"),actual=.2,forecast=.19)
  x<-list(forecasts=fc,parameters=data.frame(),cfg=c2_config("smoke"),origins=500L,
    audit=data.frame(),failures=data.frame(dimension=1L,window=500L,origin_date=as.Date("2010-01-01"),reason="test failure"))
  out<-tempfile();dir.create(out);on.exit(unlink(out,recursive=TRUE))
  expect_silent(c2_write_empirical(x,out))
  expect_true(file.exists(file.path(out,"REPORT_EN.md")))
  expect_true(any(grepl("test failure",readLines(file.path(out,"estimation_failures.csv")))))
})

test_that("nested information sets share B's largest-set eligibility calendar", {
  panel<-c2_test_panel();names(panel)[2:6]<-c("AAPL","ALD","AMGN","AXP","BA")
  panel$ALD[100]<-NA_real_
  cfg<-c2_config("smoke");cfg$window<-80L;cfg$horizons<-1L
  cfg$sets<-list("AAPL",c("AAPL","ALD"));cfg$sensitivity_windows<-integer()
  x<-c2_run_empirical(panel,cfg,c(100L,101L))
  expect_false(any(x$forecasts$origin_index==100L))
  a<-x$forecasts[x$forecasts$dimension==1L & x$forecasts$model=="mfBm",]
  b<-x$forecasts[x$forecasts$dimension==2L & x$forecasts$model=="mfBm",]
  expect_equal(c2_forecast_key(a),c2_forecast_key(b))
  expect_true(any(x$audit$reason=="asset_missing_at_origin"))
  legacy<-x;legacy$protocol_version<-NULL
  # This additional row represents a prior per-set forecast excluded by B's rule.
  extra<-a[1,,drop=FALSE];extra$origin_index<-100L;extra$origin_date<-panel$date[100]
  legacy$forecasts<-rbind(legacy$forecasts,extra)
  aligned<-c2_apply_joint_calendar(legacy,panel)
  expect_equal(aligned$forecasts,x$forecasts,tolerance=0)
  expect_equal(aligned$calendar_alignment$removed_forecast_rows,1L)
  expect_identical(c2_apply_joint_calendar(aligned,panel),aligned)
})

test_that("C baseline agrees with B's public multi-model runner", {
  panel<-c2_test_panel();names(panel)[2:6]<-c("AAPL","ALD","AMGN","AXP","BA")
  panel$ALD[40]<-NA_real_;panel$AAPL[105]<-NA_real_
  cfg<-c2_config("smoke");cfg$window<-80L;cfg$horizons<-c(1L,5L)
  cfg$sets<-list("AAPL",c("AAPL","ALD"));cfg$sensitivity_windows<-integer()
  own<-c2_run_empirical(panel,cfg,100L)$forecasts
  bcfg<-code_b_config("full");bcfg$window<-cfg$window;bcfg$horizons<-cfg$horizons
  bcfg$origin_dates<-panel$date[100]
  # B's logger is redirected by cwd to a disposable directory, not its snapshot.
  sandbox<-tempfile();dir.create(sandbox);previous<-setwd(sandbox)
  on.exit({setwd(previous);unlink(sandbox,recursive=TRUE)})
  b<-run_all_models(panel,c("AAPL","ALD"),bcfg,model_names=c("fBm","bfBm"))
  for(d in 1:2) {
    a<-own[own$model=="mfBm" & own$dimension==d,]
    z<-b[b$model==c("fBm","bfBm")[d],]
    expect_equal(c2_forecast_key(a),c2_forecast_key(z))
    expect_equal(a$forecast,z$forecast,tolerance=1e-12)
    expect_equal(a$actual,z$actual,tolerance=0)
  }
})

test_that("a refreshed failed run cannot retain previous scientific figures", {
  fc<-data.frame(model="mfBm",dimension=1L,window=500L,h=1L,
    origin_date=as.Date("2010-01-01"),target_date=as.Date("2010-01-02"),actual=.2,forecast=.19)
  x<-list(forecasts=fc,parameters=data.frame(),cfg=c2_config("smoke"),origins=500L,
    audit=data.frame(),failures=data.frame())
  out<-tempfile();dir.create(file.path(out,"figures"),recursive=TRUE)
  on.exit(unlink(out,recursive=TRUE))
  old<-file.path(out,"figures",c("cumulative_loss.png","rolling_parameters.png"))
  invisible(file.create(old))
  c2_write_empirical(x,out)
  expect_false(any(file.exists(old)))
  expect_false(any(grepl("!\\[",readLines(file.path(out,"REPORT_EN.md")))))
})
