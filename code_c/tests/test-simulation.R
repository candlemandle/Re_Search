test_that("simulation keeps oracle applicability and exact paired outcomes visible", {
  expect_true(exists("c2_run_simulation",mode="function"))
  if(!exists("c2_run_simulation",mode="function")) return(invisible(NULL))
  cfg <- c2_config("smoke"); cfg$sim_n <- 80L; cfg$sim_reps <- 2L
  cfg$horizons <- c(1L,5L); cfg$quadrature <- 64L; cfg$fit_quadrature <- 32L
  cache <- tempfile("c2_sim_cache_"); on.exit(unlink(cache,recursive=TRUE))
  set.seed(711); before <- .Random.seed
  first <- c2_run_simulation(cfg,cache_dir=cache)
  expect_identical(.Random.seed,before)
  expect_equal(nrow(first$forecasts),3L*2L*4L*2L)
  expect_setequal(unique(first$forecasts$dgp),c("mfBm","mfOU_slow","mfOU_fast"))
  expect_named(first$failures,c("dgp","replicate","model","reason"))
  f <- first$forecasts
  for(key in unique(paste(f$dgp,f$replicate,f$h))) {
    z <- f[paste(f$dgp,f$replicate,f$h)==key,]
    expect_equal(length(unique(z$actual)),1L)
    expect_equal(length(unique(z$actual_log)),1L)
  }
  na_oracle <- f[f$dgp=="mfBm" & f$model=="mfOU_known",]
  expect_true(all(na_oracle$status=="not_applicable"))
  expect_true(all(is.na(na_oracle$forecast)))
  expect_true(all(grepl("no positive true kappa",na_oracle$reason)))
  expect_true(all(f$parameter_knowledge[grepl("estimated",f$model)]=="estimated_training_only"))
  expect_true(all(f$model_specification[f$dgp!="mfBm" & f$model=="mfBm_known_driver"]=="misspecified_covariance_known_driver"))
  p <- first$parameter_errors
  expect_true(all(grepl("estimated",p$model)))
  expect_true(all(is.na(p$true_value[p$dgp=="mfBm" & p$parameter %in% c("kappa","mu")])))
  expect_true(all(first$metrics$n_requested==2L))
  second <- c2_run_simulation(cfg,cache_dir=cache)
  expect_equal(second$n_cache_hits,6L)
  expect_identical(first$forecasts,second$forecasts)
  expect_identical(first$parameter_errors,second$parameter_errors)
  expect_identical(first$comparisons,second$comparisons)
  expect_identical(first$signature,second$signature)
  checkpoint <- list.files(cache,pattern="[.]rds$",recursive=TRUE,full.names=TRUE)[1]
  unlink(checkpoint)
  writeLines("unfinished write",paste0(checkpoint,".partial"))
  resumed <- c2_run_simulation(cfg,cache_dir=cache)
  expect_equal(resumed$n_cache_hits,5L)
  expect_identical(first$forecasts,resumed$forecasts)
  expect_identical(first$parameter_errors,resumed$parameter_errors)
})

test_that("Monte Carlo errors use replication-level squared losses and paired keys", {
  expect_true(exists("c2_simulation_metrics",mode="function"))
  if(!exists("c2_simulation_metrics",mode="function")) return(invisible(NULL))
  f <- data.frame(dgp="mfOU_fast",replicate=rep(1:2,2),h=1L,
    model=rep(c("mfBm_estimated","mfOU_estimated"),each=2),status="ok",
    actual=c(2,4,2,4),forecast=c(1,1,2,2),actual_log=log(c(2,4,2,4)),
    forecast_log=log(c(1,1,2,2)),error=c(-1,-3,0,-2),
    log_error=c(-.5,-1.5,0,-1))
  out <- c2_simulation_metrics(f,requested_reps=3L)
  m <- out$metrics[out$metrics$model=="mfBm_estimated",]
  expect_equal(m$MSFE,5)
  expect_equal(m$RMSFE,sqrt(5))
  expect_equal(m$MCSE_MSFE,4)
  expect_equal(m$MSFE_log,1.25)
  expect_equal(m$MCSE_MSFE_log,1)
  expect_equal(m$n_requested,3L); expect_equal(m$n_success,2L)
  c <- out$comparisons[out$comparisons$comparison=="estimated_model_comparison" &
    out$comparisons$loss=="MSFE",]
  expect_equal(c$n_paired,2L)
  expect_equal(c$mean_loss_difference,-3)
  expect_equal(c$MCSE_loss_difference,2)
  expect_equal(c$improvement_pct,60)
})

test_that("source or configuration changes invalidate simulation checkpoint identity", {
  expect_true(exists("c2_simulation_signature",mode="function"))
  if(!exists("c2_simulation_signature",mode="function")) return(invisible(NULL))
  cfg <- c2_config("smoke")
  a <- c2_simulation_signature(cfg,source_identity=c(a="123",b="456"))
  changed <- cfg; changed$sim_delta <- cfg$sim_delta*2
  expect_false(identical(a,c2_simulation_signature(changed,source_identity=c(a="123",b="456"))))
  expect_false(identical(a,c2_simulation_signature(cfg,source_identity=c(a="123",b="457"))))
})

test_that("conditioning failures retain fitted diagnostics without replacing forecasts", {
  cfg <- c2_config("smoke");cfg$sim_n <- 80L;cfg$sim_reps <- 1L
  cfg$horizons <- c(1L,5L);cfg$quadrature <- 64L;cfg$fit_quadrature <- 32L
  # Fault injection is isolated to a cloned function environment: no production binding changes.
  run <- c2_run_simulation
  isolated <- new.env(parent=environment(run));environment(run) <- isolated
  isolated$c2_mfou_condition <- function(...) stop("injected conditioning failure")
  out <- run(cfg)
  f <- out$forecasts[out$forecasts$model=="mfOU_estimated",]
  expect_true(all(f$status=="failed"))
  expect_true(all(is.na(f$forecast)))
  expect_true(all(grepl("injected conditioning failure",f$reason)))
  expect_true(any(out$parameter_errors$model=="mfOU_estimated"))
  expect_true(any(out$diagnostics$model=="mfOU_estimated"))
  m <- out$metrics[out$metrics$model=="mfOU_estimated",]
  expect_true(all(m$n_success==0L));expect_true(all(m$n_failed==1L))
})

test_that("failed replications and undefined parameter truth remain visible in summaries", {
  f <- data.frame(dgp="mfBm",replicate=rep(1:2,2),h=1L,
    model=rep(c("mfBm_estimated","mfOU_estimated"),each=2),
    status=c("ok","ok","ok","failed"),actual=c(2,4,2,4),
    forecast=c(1,1,2,NA_real_),actual_log=log(c(2,4,2,4)),
    forecast_log=log(c(1,1,2,NA_real_)),error=c(-1,-3,0,NA_real_),
    log_error=c(-.5,-1.5,0,NA_real_))
  p <- data.frame(dgp="mfBm",model="mfOU_estimated",replicate=1:2,asset="asset1",
    parameter="kappa",estimate=c(.2,.4),true_value=NA_real_,error=NA_real_,
    truth_scope="not_defined_for_mfBm")
  out <- c2_simulation_metrics(f,p,requested_reps=2L)
  m <- out$metrics[out$metrics$model=="mfOU_estimated",]
  expect_equal(m$n_success,1L);expect_equal(m$n_failed,1L)
  expect_true(is.na(m$MCSE_MSFE))
  x <- out$comparisons[out$comparisons$comparison=="estimated_model_comparison",]
  expect_true(all(x$n_paired==1L));expect_true(all(is.na(x$MCSE_loss_difference)))
  expect_equal(out$parameter_summary$n_with_truth,0L)
  expect_true(is.na(out$parameter_summary$bias))
  expect_true(is.na(out$parameter_summary$RMSE))
  duplicated <- rbind(f,f[1,,drop=FALSE])
  expect_error(c2_simulation_metrics(duplicated),"duplicate")
})
