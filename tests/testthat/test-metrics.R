# tests for code B losses and the model confidence set

test_that("squared error and QLIKE match hand calculations", {
  expect_equal(loss_se(c(1, 2, 3), c(2, 2, 1)), c(1, 0, 4))
  # QLIKE(1, 2) = 1/2 - log(1/2) - 1
  expect_equal(loss_qlike(1, 2), 0.5 + log(2) - 1)
  expect_equal(loss_qlike(c(0.3, 1.7), c(0.3, 1.7)), c(0, 0))
  # QLIKE is never negative and smallest at the true value
  f <- seq(0.5, 2, by = 0.1)
  expect_true(all(loss_qlike(1, f) >= 0))
  expect_equal(f[which.min(loss_qlike(1, f))], 1)
})

test_that("MSFE x 100 of evaluate_class equals the hand calculation", {
  fc <- data.frame(
    model = rep(c("fBm", "bfBm"), each = 3), h = 1,
    target_date = rep(as.Date("2010-01-04") + 0:2, 2),
    forecast = c(0.2, 0.3, 0.4, 0.25, 0.3, 0.35),
    actual = rep(c(0.25, 0.32, 0.38), 2)
  )
  tab <- evaluate_class(fc, c("fBm", "bfBm"), c("2010-01-01", "2010-12-31"), "MSFE", B = 50, block = 2)
  by_hand <- c(mean(c(0.05, 0.02, 0.02)^2), mean(c(0, 0.02, 0.03)^2)) * 100
  expect_equal(tab$value, by_hand)
  expect_equal(tab$n_forecasts, c(3, 3))
  # the period filter uses the target date
  tab2 <- evaluate_class(fc, c("fBm", "bfBm"), c("2010-01-05", "2010-12-31"), "MSFE", B = 50, block = 2)
  expect_equal(tab2$n_forecasts, c(2, 2))
})

test_that("block bootstrap indices are valid and come in consecutive blocks", {
  set.seed(1)
  idx <- block_bootstrap_index(103, 10)
  expect_length(idx, 103)
  expect_true(all(idx >= 1 & idx <= 103))
  # inside a block the next index is +1 (or wraps from 103 to 1)
  step <- diff(idx[1:10])
  expect_true(all(step == 1 | step == -102))
})

test_that("MCS keeps the best model and drops a clearly worse one", {
  set.seed(2)
  n <- 500
  good <- rchisq(n, 1)
  L <- cbind(good = good, same = good + rnorm(n, 0, 0.01), bad = good + 1)
  p <- mcs_pvalues(L, B = 300, block = 5)
  expect_true(all(p >= 0 & p <= 1))
  expect_lt(p[3], 0.01) # the bad model is out
  expect_gt(max(p[1:2]), 0.5) # the best one has a large p-value
  expect_equal(max(p), 1)
})

test_that("MCS p-values are reproducible and identical models do not break it", {
  set.seed(3)
  L <- matrix(rchisq(300 * 3, 1), 300, 3)
  L[, 2] <- L[, 1] # two identical models
  set.seed(9); p1 <- mcs_pvalues(L, B = 200, block = 10)
  set.seed(9); p2 <- mcs_pvalues(L, B = 200, block = 10)
  expect_equal(p1, p2)
  expect_false(anyNA(p1))
  expect_error(mcs_pvalues(cbind(1:10, c(NA, 2:10))), "finite")
})

test_that("MCS is reported as unavailable when the block covers the whole sample", {
  # with block >= n every circular resample is a rotation of the full sample, so all
  # bootstrap means equal the sample mean; p-values of 1 would carry no information
  set.seed(4)
  L <- cbind(a = rchisq(9, 1), b = rchisq(9, 1) + 5)
  idx <- block_bootstrap_index(9, 20)
  expect_equal(sort(idx), 1:9)
  expect_true(all(is.na(mcs_pvalues(L, B = 50, block = 20))))
  expect_equal(mcs_status(9, 20), "unavailable_block_ge_n")
  expect_equal(mcs_status(1, 1), "unavailable_n_lt_2")
  expect_equal(mcs_status(21, 20), "ok")
  fc <- data.frame(model = rep(c("a", "b"), each = 9), h = 1,
                   target_date = rep(as.Date("2010-01-04") + 0:8, 2),
                   forecast = 1, actual = c(L[, 1], L[, 2]))
  tab <- evaluate_class(fc, c("a", "b"), c("2010-01-01", "2010-12-31"), "MSFE", B = 50, block = 20)
  expect_true(all(is.na(tab$mcs_p)))
  expect_equal(unique(tab$mcs_status), "unavailable_block_ge_n")
  expect_equal(unique(tab$mcs_block), 20)
})
