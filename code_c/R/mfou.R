# Stationary first-kind mfOU, eta=0 driver. No Euler/AR approximation.
# Cov(X_i(t),X_j(s)) = C_ij(t-s); unequal kappas can break reversibility.
c2_gauss_rule <- local({
  cache <- new.env(parent=emptyenv())
  function(n=128L) {
    key <- as.character(n)
    if (!exists(key,cache,inherits=FALSE)) {
      j <- seq_len(n-1L); b <- j/sqrt(4*j*j-1)
      J <- matrix(0,n,n); J[cbind(j,j+1L)] <- b; J[cbind(j+1L,j)] <- b
      e <- eigen(J,symmetric=TRUE)
      assign(key,list(z=(e$values+1)/2,w=e$vectors[1,]^2),cache)
    }
    get(key,cache,inherits=FALSE)
  }
})

c2_mfou_params <- function(H,kappa,sigma,rho,mu=rep(0,length(H))) {
  d <- length(H); rho <- as.matrix(rho)
  if (any(!is.finite(c(H,kappa,sigma,mu))) || any(H<=0|H>=1) ||
      length(kappa)!=d || length(sigma)!=d || length(mu)!=d ||
      any(kappa<=0|sigma<=0) || !identical(dim(rho),c(d,d)) ||
      any(!is.finite(rho)) || max(abs(rho-t(rho)))>1e-10 || any(abs(diag(rho)-1)>1e-10))
    stop("invalid mfOU parameters")
  noise <- mfbm_params(H,rho,0,sigma)
  if (!is_admissible_mfbm(noise)) stop("mfOU driver parameters are not admissible")
  list(H=H,kappa=kappa,sigma=sigma,rho=rho,mu=mu,d=d)
}

c2_mfou_cross_cov <- function(tau,par,i,j,quadrature=128L) {
  tau <- as.numeric(tau)
  if (any(!is.finite(tau))) stop("nonfinite covariance lag")
  a <- par$H[i]+par$H[j]; ki <- par$kappa[i]; kj <- par$kappa[j]
  factor <- par$rho[i,j]*par$sigma[i]*par$sigma[j]/2
  if (factor==0) return(rep(0,length(tau)))
  out <- numeric(length(tau)); neg <- tau<0
  if (any(neg)) out[neg] <- c2_mfou_cross_cov(-tau[neg],par,j,i,quadrature)
  zero <- tau==0
  out[zero] <- factor*gamma(a+1)*(ki^(1-a)+kj^(1-a))/(ki+kj)
  pos <- tau>0
  if (any(pos)) {
    t <- tau[pos]
    # Exact Brownian anchor avoids cancellation at very large lags.
    if (abs(a-1)<1e-12) {
      out[pos] <- 2*factor*exp(-ki*t)/(ki+kj)
    } else {
      xj <- kj*t; xi <- ki*t
      plus <- a*kj^(-a)*exp(xj+lgamma(a)+pgamma(xj,a,lower.tail=FALSE,log.p=TRUE))
      rule <- c2_gauss_rule(quadrature)
      integral <- drop(exp(-outer(xi,1-rule$z^(1/a))) %*% rule$w)
      minus <- gamma(a+1)*ki^(-a)*exp(-xi)-t^a*integral
      out[pos] <- factor*(kj*plus+ki*minus)/(ki+kj)
    }
  }
  out
}

c2_mfou_cov_matrix <- function(par,times,quadrature=128L) {
  n <- length(times); d <- par$d
  dt <- outer(times,times,"-"); lags <- sort(unique(as.vector(dt)))
  idx <- matrix(match(dt,lags),n,n)
  S <- matrix(0,n*d,n*d)
  for (i in seq_len(d)) for (j in seq_len(d)) {
    cv <- c2_mfou_cross_cov(lags,par,i,j,quadrature)
    S[seq(i,n*d,by=d),seq(j,n*d,by=d)] <- matrix(cv[idx],n,n)
  }
  (S+t(S))/2
}

# Expected raw variogram and variance after subtracting the sample mean.
# The latter subtracts Var(sample mean) exactly, including real trading gaps.
c2_moment_design <- function(days,lags) {
  all_lags<-seq.int(0L,max(c(lags,max(days)-min(days))))
  count<-tabulate(abs(outer(days,days,"-"))+1L,nbins=length(all_lags))
  list(all_lags=all_lags,count=count,n=length(days))
}

c2_mfou_moment_shape <- function(H,kappa,lags,days,delta,quadrature=64L,design=NULL) {
  p <- list(H=H,kappa=kappa,sigma=1,rho=matrix(1),d=1)
  v <- c2_mfou_cross_cov(0,p,1,1,quadrature)
  if(is.null(design)) design<-c2_moment_design(days,lags)
  all_lags<-design$all_lags
  cv <- c2_mfou_cross_cov(all_lags*delta,p,1,1,quadrature)
  centered <- v-sum(design$count*cv)/design$n^2
  c(2*(v-cv[match(lags,all_lags)]),centered)
}

c2_fit_mfou_moments <- function(emp,lags,days,delta,H_start=.2,
                                kappa_starts=c(.1,1,10),quadrature=64L) {
  if (any(!is.finite(emp)|emp<=0) || length(emp)!=length(lags)+1L)
    stop("invalid mfOU empirical moments")
  lower <- c(.02,log(.02)); upper <- c(.8,log(50))
  design<-c2_moment_design(days,lags)
  objective <- function(x,details=FALSE) {
    shape <- c2_mfou_moment_shape(x[1],exp(x[2]),lags,days,delta,quadrature,design)
    z <- shape/emp
    if (any(!is.finite(z)|z<=0)) return(if(details) NULL else 1e30)
    s2 <- sum(z)/sum(z^2)
    value <- mean((s2*z-1)^2)
    if(details) list(value=value,sigma2=s2) else value
  }
  fits <- lapply(kappa_starts,function(k) tryCatch(
    optim(c(min(max(H_start,.03),.79),log(k)),objective,method="L-BFGS-B",
      lower=lower,upper=upper,control=list(maxit=120, factr=1e7)),error=function(e) NULL))
  ok <- vapply(fits,function(x) !is.null(x) && x$convergence==0L && is.finite(x$value),TRUE)
  if(!any(ok)) stop("mfOU moment estimation failed at every start")
  successful <- fits[ok]; best <- successful[[which.min(vapply(successful,`[[`,0,"value"))]]
  detail <- objective(best$par,TRUE); H <- best$par[1]; k <- exp(best$par[2]); s2 <- detail$sigma2
  theta <- c(H,log(k),log(s2))
  moment <- function(v) exp(v[3])*c2_mfou_moment_shape(v[1],exp(v[2]),lags,days,delta,quadrature,design)/emp
  J <- sapply(seq_len(3),function(j) {
    up <- down <- theta; up[j] <- up[j]+1e-4; down[j] <- down[j]-1e-4
    (moment(up)-moment(down))/2e-4
  })
  sv <- svd(J)$d
  list(H=H,kappa=k,sigma2=s2,objective=detail$value,n_starts=length(fits),n_converged=sum(ok),
    boundary=any(abs(best$par-lower)<1e-4 | abs(best$par-upper)<1e-4),
    jacobian_rank=sum(sv>max(sv)*1e-6),jacobian_condition=max(sv)/min(sv),
    kappa_start_spread=max(vapply(successful,function(f) exp(f$par[2]),0))/
      min(vapply(successful,function(f) exp(f$par[2]),0)),
    kappa_window_span=k*(max(days)-min(days))*delta)
}

c2_fit_mfou <- function(Y,delta,lags=c(1L,2L,5L,10L,20L,60L),quadrature=64L) {
  Y <- as.matrix(Y);Y[!is.finite(Y)] <- NA_real_
  if(is.null(colnames(Y))) colnames(Y)<-paste0("asset",seq_len(ncol(Y)))
  obs <- which(complete.cases(Y)); days <- obs-1L
  if(length(obs)<30L) stop("too few joint observations for mfOU")
  pilot <- window_mfbm_params_calendar(Y,delta)
  fits <- lapply(seq_len(ncol(Y)),function(i) {
    emp <- vapply(lags,function(l) {
      z <- Y[(l+1):nrow(Y),i]-Y[1:(nrow(Y)-l),i]
      jointly <- complete.cases(Y[(l+1):nrow(Y),,drop=FALSE]) & complete.cases(Y[1:(nrow(Y)-l),,drop=FALSE])
      mean(z[jointly]^2)
    },0)
    emp <- c(emp,mean((Y[obs,i]-mean(Y[obs,i]))^2))
    c2_fit_mfou_moments(emp,lags,days,delta,pilot$H[i],quadrature=quadrature)
  })
  H <- vapply(fits,`[[`,0,"H"); k <- vapply(fits,`[[`,0,"kappa"); sigma <- sqrt(vapply(fits,`[[`,0,"sigma2"))
  d <- ncol(Y); rho <- diag(d); increments <- diff(Y)
  increments <- increments[complete.cases(increments),,drop=FALSE]
  for(i in seq_len(d)) for(j in seq_len(d)) if(i<j) {
    unit <- list(H=H,kappa=k,sigma=sigma,rho=matrix(1,d,d),d=d)
    cov_inc <- 2*c2_mfou_cross_cov(0,unit,i,j)-c2_mfou_cross_cov(delta,unit,i,j)-c2_mfou_cross_cov(-delta,unit,i,j)
    rho[i,j] <- rho[j,i] <- mean(increments[,i]*increments[,j])/cov_inc
  }
  admissible <- admissible_window_params(list(H=H,sigma2=sigma^2,rho=rho))
  p <- c2_mfou_params(admissible$H,k,sigma,admissible$rho,colMeans(Y[obs,,drop=FALSE]))
  diagnostics <- do.call(rbind,lapply(seq_len(d),function(i) data.frame(asset=colnames(Y)[i],
    H=H[i],kappa=k[i],sigma=sigma[i],objective=fits[[i]]$objective,
    boundary=fits[[i]]$boundary,n_starts=fits[[i]]$n_starts,n_converged=fits[[i]]$n_converged,
    jacobian_rank=fits[[i]]$jacobian_rank,jacobian_condition=fits[[i]]$jacobian_condition,
    kappa_start_spread=fits[[i]]$kappa_start_spread,kappa_window_span=fits[[i]]$kappa_window_span,
    rho_shrink_steps=admissible$shrink,drift_half_life_days=log(2)/k[i]/delta)))
  list(par=p,diagnostics=diagnostics)
}

c2_mfou_condition <- function(Y,horizons,delta,par,estimate_mean=TRUE,quadrature=128L) {
  Y <- as.matrix(Y); obs <- which(complete.cases(Y)); n <- nrow(Y); d <- ncol(Y)
  if(!length(obs) || tail(obs,1)!=n || d!=par$d) stop("invalid mfOU conditioning observations")
  times <- (obs-1)*delta; S <- c2_mfou_cov_matrix(par,times,quadrature)
  R <- chol(S) # no hidden diagonal jitter or fallback
  solveS <- function(z) backsolve(R,forwardsolve(t(R),z))
  y <- as.vector(t(Y[obs,,drop=FALSE])); D <- kronecker(matrix(1,length(obs),1),diag(d))
  mu <- par$mu
  if(estimate_mean) mu <- drop(solve(crossprod(D,solveS(D)),crossprod(D,solveS(y))))
  gamma_xy <- matrix(0,length(y),length(horizons))
  for(j in seq_len(d)) {
    lag <- outer(times,(n-1+horizons)*delta,"-")
    gamma_xy[seq(j,length(y),by=d),] <- matrix(c2_mfou_cross_cov(lag,par,j,1,quadrature),length(obs))
  }
  W <- solveS(gamma_xy)
  m <- mu[1]+drop(crossprod(y-drop(D%*%mu),W))
  v <- c2_mfou_cross_cov(0,par,1,1,quadrature)-colSums(gamma_xy*W)
  if(any(v < -1e-8) || any(!is.finite(c(m,v)))) stop("invalid mfOU conditional moments")
  v <- pmax(v,0)
  list(forecast=exp(m+v/2),log_mean=m,log_variance=v,mu=mu,
    drift_half_life_days=log(2)/par$kappa/delta,
    cholesky_diagonal_ratio=max(diag(R))/min(diag(R)),n_joint_observations=length(obs))
}
