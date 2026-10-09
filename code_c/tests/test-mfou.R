test_that("stationary covariance has the analytic multivariate OU limit", {
  p <- c2_mfou_params(c(.5,.5),c(.7,1.2),c(.8,1.1),matrix(c(1,.4,.4,1),2),c(-2,-1))
  t <- c(0,.001,.1,1,10)
  expect_equal(c2_mfou_cross_cov(t,p,1,2), .4*.8*1.1*exp(-.7*t)/1.9,tolerance=1e-9)
  expect_equal(c2_mfou_cross_cov(-t,p,1,2), .4*.8*1.1*exp(-1.2*t)/1.9,tolerance=1e-9)
})

test_that("rough covariance agrees with independent exponential integrals", {
  p <- c2_mfou_params(c(.1,.2),c(.7,1.2),c(1,.9),matrix(c(1,.4,.4,1),2),c(0,0))
  a <- .3; tau <- .4
  plus <- integrate(function(u) 1.2*exp(-1.2*u)*(tau+u)^a,0,Inf,rel.tol=1e-10)$value
  minus <- integrate(function(u) .7*exp(-.7*u)*abs(tau-u)^a,0,tau,rel.tol=1e-10)$value +
    integrate(function(u) .7*exp(-.7*u)*abs(tau-u)^a,tau,Inf,rel.tol=1e-10)$value
  expected <- .4*.9/2*(1.2/1.9*plus+.7/1.9*minus-tau^a)
  expect_equal(c2_mfou_cross_cov(tau,p,1,2),expected,tolerance=2e-7)
  expect_equal(c2_mfou_cross_cov(0,p,1,1),gamma(1.2)/(2*.7^.2),tolerance=1e-12)
  S <- c2_mfou_cov_matrix(p,seq(0,1,length.out=30))
  expect_equal(S,t(S),tolerance=1e-12)
  expect_gt(min(eigen(S,symmetric=TRUE,only.values=TRUE)$values),0)
  expect_equal(S,c2_mfou_cov_matrix(p,seq(0,1,length.out=30)+100),tolerance=1e-10)
  expect_equal(c2_mfou_cross_cov(.4,p,1,2,128),c2_mfou_cross_cov(.4,p,1,2,256),tolerance=1e-7)
})

test_that("Gaussian conditioning uses history and log-normal adjustment", {
  p <- c2_mfou_params(.5,2,.4,matrix(1),-2)
  Y <- matrix(c(-2.1,-1.9,-1.7),ncol=1)
  x <- c2_mfou_condition(Y,c(1L,5L),1/252,p,estimate_mean=FALSE)
  m <- -2+exp(-2*c(1,5)/252)*(-1.7+2)
  v <- .4^2/4*(1-exp(-4*c(1,5)/252))
  expect_equal(x$log_mean,m,tolerance=1e-9)
  expect_equal(x$log_variance,v,tolerance=1e-9)
  expect_equal(x$forecast,exp(m+v/2),tolerance=1e-9)
  expect_equal(x$drift_half_life_days,rep(log(2)/2*252,1),tolerance=1e-12)
})

test_that("population moments identify interior mfOU parameters", {
  lags <- c(1,2,5,10,20,60)
  days <- 0:249; H <- .2; k <- 3; s2 <- .7
  emp <- c2_mfou_moment_shape(H,k,lags,days,1/252)*s2
  fit <- c2_fit_mfou_moments(emp,lags,days,1/252,H_start=.25,kappa_starts=c(.1,1,10))
  expect_lt(fit$objective,1e-7)
  expect_equal(fit$H,H,tolerance=3e-3)
  expect_equal(fit$kappa,k,tolerance=.03)
  expect_equal(fit$sigma2,s2,tolerance=3e-3)
  expect_equal(fit$jacobian_rank,3L)
})

test_that("cached moment design retains gaps and exact pair multiplicities", {
  days<-c(0L,1L,3L,4L);lags<-c(1L,2L)
  design<-c2_moment_design(days,lags)
  expect_equal(sum(design$count),length(days)^2)
  expect_equal(design$count[match(0L,design$all_lags)],4L)
  expect_equal(c2_mfou_moment_shape(.2,1,lags,days,1/252,design=design),
    c2_mfou_moment_shape(.2,1,lags,days,1/252),tolerance=0)
})
