# =============================================================================
# Title: Inference for Multivariate Fractional Brownian Motion
# Source: "Modeling and Forecasting Realized Volatility with Multivariate Fractional Brownian Motion"
#         by Markus Bibinger, Jun Yu, and Chen Zhang
# Description:
# This R script implements a set of functions used to compute the standard errors 
# of key estimators in the context of modeling realized volatility using 
# multivariate fractional Brownian motion (mfBm).


## ---------------------------
## Prepared Functions
## ---------------------------

varRV<-function(n,sigm,H){
  gamma<-((1:n)+1)^(2*H)+((1:n)-1)^(2*H)-2*(1:n)^(2*H)
  varRV<-(2+sum((1-(1:n)/n)*gamma^2))*sigm^4}

varRCov<-function(n,rho,sigm1,sigm2,H1,H2){
  gamma<-((1:n)+1)^(H1+H2)+((1:n)-1)^(H1+H2)-2*(1:n)^(H1+H2)
  gamma1<-((1:(n))+1)^(2*H1)+((1:(n))-1)^(2*H1)-2*(1:(n))^(2*H1)
  gamma2<-((1:(n))+1)^(2*H2)+((1:(n))-1)^(2*H2)-2*(1:(n))^(2*H2)
  varRCov<-sigm1^2*sigm2^2*((1+rho^2)+rho^2/2*sum((1-(1:n)/n)*gamma^2)+0.5*sum((1-(1:n)/n)*gamma1*gamma2))}

covRVRCov<-function(n,rho,sigm1,sigm2,H1,H2){
  gamma<-((1:n)+1)^(H1+H2)+((1:n)-1)^(H1+H2)-2*(1:n)^(H1+H2)
  gamma1<-((1:(n))+1)^(2*H1)+((1:(n))-1)^(2*H1)-2*(1:(n))^(2*H1)
  gamma2<-((1:(n))+1)^(2*H2)+((1:(n))-1)^(2*H2)-2*(1:(n))^(2*H2)
  covRVRCov<-2*rho*sigm1^3*sigm2*(1+0.5*sum((1-(1:n)/n)*gamma1*gamma))}

covRVRV<-function(n,rho,sigm1,sigm2,H1,H2){
  gamma<-((1:n)+1)^(H1+H2)+((1:n)-1)^(H1+H2)-2*(1:n)^(H1+H2)
  covRVRV<-2*rho^2*sigm1^2*sigm2^2*(1+0.5*sum((1-(1:n)/n)*gamma*gamma))}

# Asymptotic variance of rho
varrho<-function(n,rho,sigm1,sigm2,H1,H2){
  varrho<-rho^2/(4*sigm1^4)*varRV(n,sigm1,H1)+rho^2/(4*sigm2^4)*varRV(n,sigm2,H2)+1/(sigm1^2*sigm2^2)*varRCov(n,rho,sigm1,sigm2,H1,H2)+rho^2/(2*sigm1^2*sigm2^2)*covRVRV(n,rho,sigm1,sigm2,H1,H2)-rho/(sigm1^3*sigm2)*covRVRCov(n,rho,sigm1,sigm2,H1,H2)-rho/(sigm1*sigm2^3)*covRVRCov(n,rho,sigm2,sigm1,H2,H1)}

# Simpler version of varrho (asymptotically equal)
varrho2<-function(n,rho,H1,H2){
  gamma1<-((1:(n-1))+1)^(2*H1)+((1:(n-1))-1)^(2*H1)-2*(1:(n-1))^(2*H1)
  gamma2<-((1:(n-1))+1)^(2*H2)+((1:(n-1))-1)^(2*H2)-2*(1:(n-1))^(2*H2)
  gamma12<-((1:(n-1))+1)^(H1+H2)+((1:(n-1))-1)^(H1+H2)-2*(1:(n-1))^(H1+H2)
  varrho<-(1-rho^2)^2+rho^2*(0.25*sum((1-(2:(n))/n)*(gamma1^2+gamma2^2+2*(1+rho^2)*gamma12^2)))-rho^2*sum((1-(2:(n))/n)*(gamma1*gamma12+gamma2*gamma12))+0.5*sum((1-(2:(n))/n)*gamma1*gamma2)
  print(varrho)}

# Asymptotic variance of H
varHest<-function(n,H){
  gamma<-((1:n)+1)^(2*H)+((1:n)-1)^(2*H)-2*(1:n)^(2*H)
  varHest<-1/4*1/(log(2)^2)*(4+sum((1-(1:n)/n)*gamma^2)-2*(2^(2*H)+2^(-2*H)*sum((1-(2:(n-1))/n)*((3:(n))^(2*H)-(1:(n-2))^(2*H)-(2:(n-1))^(2*H)+(0:(n-3))^(2*H))^2))+2^(-4*H)*sum((1-(1:(n-2))/n)*(-2*(1:(n-2))^(2*H)+(3:n)^(2*H)+(abs((-1):(n-4)))^(2*H))^2))
}

# Asymptotic variance of sigma^2
varsigma<-function(n,sigm,H){
  varsigma<-varHest(n,H)*4*sigm^4}

# Asymptotic variance of eta
vareta<-function(n,rho,H1,H2){
  gamma1<-((1:(n-1))+1)^(2*H1)+((1:(n-1))-1)^(2*H1)-2*(1:(n-1))^(2*H1)
  gammas1<-((2:(n))+1)^(2*H1)+((2:(n))-1)^(2*H1)-2*(2:(n))^(2*H1)
  gammass1<-c(2,((1:(n-2))+1)^(2*H1)+((1:(n-2))-1)^(2*H1)-2*(1:(n-2))^(2*H1))
  gamma2<-((1:(n-1))+1)^(2*H2)+((1:(n-1))-1)^(2*H2)-2*(1:(n-1))^(2*H2)
  gammas2<-((2:(n))+1)^(2*H2)+((2:(n))-1)^(2*H2)-2*(2:(n))^(2*H2)
  gammass2<-c(2,((1:(n-2))+1)^(2*H2)+((1:(n-2))-1)^(2*H2)-2*(1:(n-2))^(2*H2))
  gamma12<-((1:(n-1))+1)^(H1+H2)+((1:(n-1))-1)^(H1+H2)-2*(1:(n-1))^(H1+H2)
  gammas12<-((2:(n))+1)^(H1+H2)+((2:(n))-1)^(H1+H2)-2*(2:(n))^(H1+H2)
  gammass12<-c(2,((1:(n-2))+1)^(H1+H2)+((1:(n-2))-1)^(H1+H2)-2*(1:(n-2))^(H1+H2))
  avar<-2*(1-rho^2)-0.5*(2^(2*H1)-2)*(2^(2*H2)-2)+0.5*rho^2*(2^(H1+H2)-2)^2+rho^2*sum((1-(2:(n))/n)*(gammas12*gammass12-gamma12^2))+0.5*sum((1-(2:(n))/n)*(2*gamma1*gamma2-gammas1*gammass2-gammass1*gammas2))
  vareta<-avar/(2^(H1+H2)-2)^2
  print(vareta)}


## ---------------------------
## Inference
## ---------------------------

# parameters
H1<-0.1
H2<-0.4
rho<-0.4

# Truncation
n<-2000

# Sample size and sampling interval 
N<-1000
delta<-1/250

# Initialize storage
seH1     <- numeric(1)
seH2     <- numeric(1)
sesigma1 <- numeric(1)
sesigma2 <- numeric(1)
serho    <- numeric(1)
seeta    <- numeric(1)

# Standard error
seH1[1]<-(varHest(n,H1)/N)^.5
seH2[1]<-(varHest(n,H2)/N)^.5
sesigma1[1] <- (varsigma(n,1,H1)/N*log(1/delta)*log(1/delta))^.5
sesigma2[1] <- (varsigma(n,1,H2)/N*log(1/delta)*log(1/delta))^.5
serho[1] <-(varrho2(n,rho,H1,H2)/N)^.5
seeta[1] <-(vareta(n,rho,H1,H2)/N)^.5