# Run the authors' R files (author_code/BYZ, from
# https://fba.um.edu.mo/wp-content/uploads/2025/04/BYZ.zip) and compare them
# with Code A on identical inputs.
# Usage: Rscript scripts/00_check_author_code.R
# Output: results/tables/code_a_author_code_comparison.csv

source("R/code_a_config.R")
source_code_a()
cfg <- code_a_config("full")
code_a_log("00_check_author_code: start")
auth_dir <- file.path("author_code", "BYZ")
stopifnot(file.exists(file.path(auth_dir, "Setup code.R")))

# The author files run top-level code; evaluate each in its own environment.
auth <- new.env()
sys.source(file.path(auth_dir, "Setup code.R"), envir = auth)
set.seed(cfg$seed)
est_env <- new.env(parent = auth)
utils::capture.output(sys.source(file.path(auth_dir, "BYZestimator.R"), envir = est_env))
inf_env <- new.env(parent = auth)
utils::capture.output(sys.source(file.path(auth_dir, "BYZinference.R"), envir = inf_env))
code_a_log("  author files ran: BYZestimator.R (H1hat = ", round(est_env$H1hat, 4),
           "), BYZinference.R (seeta = ", round(inf_env$seeta, 4), ")")

rows <- list()
add <- function(check, quantity, author, ours, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(check = check, quantity = quantity, author = author,
                                         ours = ours, abs_diff = abs(author - ours), note = note)
}

# 1. Standard errors (BYZinference.R functions, its truncation n = 2000, N = 1000, Delta = 1/250)
q <- function(expr) { utils::capture.output(v <- expr); v } # author functions print()
for (rho in c(0, 0.4)) for (N in c(500, 1000)) {
  add("SE (BYZinference.R)", sprintf("se H1, N=%d", N), sqrt(inf_env$varHest(2000, 0.1) / N),
      sqrt(avar_hurst(0.1) / N), "author: finite truncation with (1 - r/n) weights")
  add("SE (BYZinference.R)", sprintf("se rho, rho=%.1f N=%d", rho, N),
      sqrt(q(inf_env$varrho2(2000, rho, 0.1, 0.4)) / N), sqrt(avar_rho(0.1, 0.4, rho) / N))
  add("SE (BYZinference.R)", sprintf("se eta, rho=%.1f N=%d", rho, N),
      sqrt(q(inf_env$vareta(2000, rho, 0.1, 0.4)) / N), sqrt(avar_eta(0.1, 0.4, rho) / N))
}
add("SE (BYZinference.R)", "se sigma2_1, N=1000, Delta=1/250", inf_env$sesigma1,
    se_sigma2(0.1, 1, 1000, 1 / 250))

# 2. Point estimators on the same path (author code as in BYZestimator.R)
author_byz <- function(b1, b2, delta) {
  n <- length(b1)
  H1 <- 0.5 / log(2) * log(sum((b1[3:n] - b1[1:(n - 2)])^2) / sum(diff(b1)^2))
  H2 <- 0.5 / log(2) * log(sum((b2[3:n] - b2[1:(n - 2)])^2) / sum(diff(b2)^2))
  eta <- mean(diff(b2)[2:(n - 1)] * diff(b1)[1:(n - 2)] - diff(b1)[2:(n - 1)] * diff(b2)[1:(n - 2)]) /
    (sqrt(mean(diff(b1, lag = 2)^2) * mean(diff(b2, lag = 2)^2)) - 2 * sqrt(mean(diff(b1)^2) * mean(diff(b2)^2)))
  rho <- mean(diff(b1) * diff(b2)) / sqrt(mean(diff(b1)^2) * mean(diff(b2)^2))
  c(H1 = H1, H2 = H2, s1 = mean(diff(b1)^2) / delta^(2 * H1), rho = rho, eta = eta)
}
set.seed(cfg$seed + 1)
x <- simulate_mfbm(1000, c(0.1, 0.4), 0.4, 0.3, delta = 1 / 250)[[1]]
B <- increments_to_levels(x)
a <- author_byz(B[, 1], B[, 2], 1 / 250)
e <- estimate_mfbm(x, 1 / 250)
add("BYZ estimator (same path)", "H1", a["H1"], e$H[1])
add("BYZ estimator (same path)", "H2", a["H2"], e$H[2])
add("BYZ estimator (same path)", "sigma2_1", a["s1"], e$sigma2[1],
    "author mean() uses n-1 sums: identical up to the (n-1)/n convention")
add("BYZ estimator (same path)", "rho", a["rho"], e$rho[1, 2])
add("BYZ estimator (same path)", "eta: author (8b) vs ours, printed convention", a["eta"],
    eta_mm(x[, 1], x[, 2], "printed"), "mean() normalisations differ by (n-2)/(n-1) factors")
add("BYZ estimator (same path)", "eta: author (8b) vs ours, proof convention", a["eta"],
    eta_mm(x[, 1], x[, 2]), "opposite sign: conventions differ, see 3.")

# 3. Sign convention of eta: simulate with the authors' simMFBM, eta = 0.5
set.seed(cfg$seed + 2)
sims <- suppressWarnings({
  utils::capture.output(z <- auth$simMFBM(n = 2000, H = c(0.1, 0.4), sig = c(1, 1),
                                         rho = matrix(c(1, .4, .4, 1), 2),
                                         eta = matrix(c(0, .5, -.5, 0), 2), print = FALSE, nbSamp = 50))
  z
})
eta_author_sim <- vapply(sims, function(s) eta_mm(diff(s[, 1]), diff(s[, 2]), "printed"), 0)
eta_ours_sim <- vapply(sims, function(s) eta_mm(diff(s[, 1]), diff(s[, 2]), "proof"), 0)
add("eta convention (simMFBM, eta = 0.5, 50 paths)", "mean of (8b) as printed / author code",
    0.5, mean(eta_author_sim), "author simulator + printed (8b) recover +eta")
add("eta convention (simMFBM, eta = 0.5, 50 paths)", "mean of Appendix C.2 version",
    0.5, mean(eta_ours_sim), "= -eta: simMFBM uses (rho - eta sign) i.e. eta of opposite sign to eq. (4)")
set.seed(cfg$seed + 3)
ours <- simulate_mfbm(2000, c(0.1, 0.4), 0.4, 0.5, delta = 1, nsim = 50)
add("eta convention (our simulator, eta = 0.5 in eq. (4), 50 paths)", "mean of Appendix C.2 version",
    0.5, mean(vapply(ours, function(s) eta_mm(s[, 1], s[, 2]), 0)))

# 4. Amblard-Coeurjolly estimator: author estMFBM vs our estimate_ac on the same path
set.seed(cfg$seed + 4)
x <- simulate_mfbm(1000, c(0.1, 0.4), 0.4, 0.5, delta = 1 / 250)[[1]]
B <- increments_to_levels(x) # same levels (incl. B_0 = 0) for both implementations
ac_a <- auth$estMFBM(B, nma = "i2", M1 = 1, M2 = 5, w = c(1, 0, 0), forceH = FALSE, naiveRhoEta = FALSE)
ac_o <- estimate_ac(x, 1 / 250)
add("AC estimator (same path)", "H1", ac_a$H[1], ac_o$H[1])
add("AC estimator (same path)", "H2", ac_a$H[2], ac_o$H[2])
add("AC estimator (same path)", "rho", ac_a$rho[1, 2], ac_o$rho[1, 2])
add("AC estimator (same path)", "|eta|", ac_a$eta[1, 2], ac_o$eta[1, 2])

out <- do.call(rbind, rows)
rownames(out) <- NULL
code_a_write(out, "code_a_author_code_comparison.csv")
print(out[, c("check", "quantity", "author", "ours", "abs_diff")], digits = 4)
code_a_log_run("00_check_author_code", "full", cfg)
code_a_log("00_check_author_code: done")
