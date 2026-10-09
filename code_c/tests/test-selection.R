test_that("asset selection uses reviewed calendar predictions without future inputs", {
  set.seed(16); n <- 130
  panel <- data.frame(date=as.Date("2010-01-01")+1:n,
    matrix(cumsum(rnorm(n*5,sd=.03)),n,5))
  names(panel)[2:6] <- c("AAPL","ALD","AMGN","AXP","BA")
  cfg <- c2_config("smoke"); cfg$window <- 80L; cfg$horizons <- c(1L,5L)
  cfg$selection_targets <- "AAPL";cfg$selection_universe <- names(panel)[-1]
  cfg$selection_k <- 1L;cfg$selection_fixed <- list(AAPL="ALD")
  panel$ALD[35] <- NA_real_
  x <- c2_run_selection(panel,cfg,100L)
  f <- x$forecasts[x$forecasts$strategy=="fixed",]
  expect_equal(f$forecast,forecast_mfbm_calendar(as.matrix(panel[21:100,c("AAPL","ALD")]),c(1L,5L),1/252),tolerance=1e-12)
  z <- panel;z$ALD[101:n] <- NA_real_
  y <- c2_run_selection(z,cfg,100L)
  expect_equal(x$forecasts,y$forecasts,tolerance=1e-12)
  expect_equal(x$selected_assets,y$selected_assets,tolerance=1e-12)
  expect_true(all(x$forecasts$target_date==panel$date[100L+x$forecasts$h]))
})

test_that("window sensitivity compares losses only on intersected forecast keys", {
  fc <- expand.grid(model=c("mfBm","mfOU"),window=c(250L,500L,1000L),day=1:4)
  fc$dimension <- 1L;fc$h <- 1L;fc$origin_date <- as.Date("2020-01-01")+fc$day
  fc$target_date <- fc$origin_date+1;fc$actual <- .2;fc$forecast <- .19
  fc <- fc[!(fc$window==1000L & fc$day==1),]
  x <- c2_window_comparison(fc)
  expect_true(all(x$n[x$alternative_window==1000L]==3L))
  expect_true(all(x$n[x$alternative_window==250L]==4L))
  expect_true(all(x$gain_pct==0))
})
