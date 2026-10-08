# =============================================================================
# Title: Estimation for Multivariate Fractional Brownian Motion
# Source: "Modeling and Forecasting Realized Volatility with Multivariate Fractional Brownian Motion"
#         by Markus Bibinger, Jun Yu, and Chen Zhang
# Description:
# This R script implements a set of functions used to obtain the estimators 
# in the context of modeling realized volatility using 
# multivariate fractional Brownian motion (mfBm).
# Before running this R code, make sure to execute the preliminary setup code.


## Simulation for Bivariate fBm Estimation

# --- Parameters ---
H1 <- 0.1
H2 <- 0.4
rhov <- 0
etav <- 0

# --- Sample size ---
n <- 500
delta <- 1  # Sampling interval

# --- Initialize storage for estimators ---
etahat     <- numeric(1)
rhohat     <- numeric(1)
H1hat      <- numeric(1)
H2hat      <- numeric(1)
sigma1hat  <- numeric(1)
sigma2hat  <- numeric(1)

# --- Generate bivariate fBm with delta = 1 ---
Z <- simMFBM(n = 500 ,H = c(H1, H2),sig = c(1, 1), 
             rho = matrix(c(1, rhov, rhov, 1), nc = 2),
             eta = matrix(c(0, etav, -etav, 0), nc = 2), 
             plot = FALSE, print = TRUE, choix = NULL,
             forceEta = FALSE,nbSamp = 1)

# --- Rescale process according to delta ---
b1 <- Z[, 1] * delta^H1
b2 <- Z[, 2] * delta^H2

# --- MYZ Estimator ---
# Equation 6: Estimate H1 and H2
H1hat[1] <- 0.5 / log(2) * log(
  sum((b1[3:n] - b1[1:(n - 2)])^2) / sum(diff(b1)^2)
)
H2hat[1] <- 0.5 / log(2) * log(
  sum((b2[3:n] - b2[1:(n - 2)])^2) / sum(diff(b2)^2)
)

# Equation 7: Estimate sigma1 and sigma2
sigma1hat[1] <- mean(diff(b1)^2) / delta^(2 * H1hat[1])
sigma2hat[1] <- mean(diff(b2)^2) / delta^(2 * H2hat[1])

# Equation 8: Estimate eta and rho
etahat[1] <- mean(
  diff(b2)[2:(n - 1)] * diff(b1)[1:(n - 2)] -
    diff(b1)[2:(n - 1)] * diff(b2)[1:(n - 2)]
) / (
  sqrt(mean(diff(b1, lag = 2)^2) * mean(diff(b2, lag = 2)^2)) -
    2 * sqrt(mean(diff(b1)^2) * mean(diff(b2)^2))
)

rhohat[1] <- mean(diff(b1) * diff(b2)) / sqrt(
  mean(diff(b1)^2) * mean(diff(b2)^2)
)

# --- AC Estimator ---
ZZ <- cbind(b1, b2)  # Combine into matrix
estMFBM(
  ZZ,
  nma = "i2",
  M1 = 1,
  M2 = 5,
  w = c(1, 0, 0),
  forceH = FALSE,
  naiveRhoEta = FALSE
)
