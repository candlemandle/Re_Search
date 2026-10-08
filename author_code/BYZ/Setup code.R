existMFBB<-function(H,rhomat,etamat=NULL){
  # Check if process (mfBB) exists for given parameters
  # H: Hurst parameter vector
  # rhomat: matrix  with pairwise correlations (at time 1)
  # etamat: matrix with asymmetry parameters
  # Output: negative eigenvalues, mfBB exists iff Output = 0
  # Based on: https://sites.google.com/site/homepagejfc/software
  p<-length(H)
  C<-matrix(0,nr=p,nc=p)
  for (j in (1:p)){for (k in (1:p)){
    Hjk<-H[j]+H[k]
    if (Hjk!=1){
      pre<-rhomat[j,k]*sin(pi/2*Hjk) 
      pim<- etamat[j,k]*cos(pi/2*Hjk)
    }
    else{
      pre<-rhomat[j,k];pim<- pi/2*etamat[j,k]
    }
    C[j,k]<-complex(real=pre,imaginary=pim)
    C[j,k]<-gamma(H[j]+H[k]+1)*C[j,k]
  }}
  ltmp<-sum(eigen(C)$values<0)
  ltmp
}

simMFBM<- function(n=50, H=c(.2,.2),
                   sig=c(1,1),rho=matrix(c(1,.5,.5,1),nc=2),eta=matrix(c(0,.5,-.5,0),nc=2),
                   plot=FALSE,print=TRUE,choix=NULL,forceEta=FALSE,nbSamp=1){ 
  
  ## Simulation of a multivariate fractional Brownian motion.
  ## nbSamp : number of sample paths if =1 a (n,p) matrix is returned
  ## if >1 a list of (n,p) matrices is returned 
  ##  notation: causal, wb, eta stand for the causal MFBM, the well-balanced MFBM, 
  ## a general MFBM with prescribed eta matrix. 
  ## forceEta=TRUE means that eta_jk is changed by eta_jk/(1-H[j]-H[k]) only when H[j]+H[k] !=1
  ## so that this is easy to check the semidefinite positiveness using ellipsis in 
  ## See the reference for details.
  ##
  ## Reference. Amblard, Coeurjolly, Philippe and Lavancier (2012)
  ## Basic properties of the multivariate fractional Brownian motion.
  ## Bulletin de la SMF, Séminaires et Congrés, 28, 65-87.
  
  p<-length(H)
  G<-matrix(0,nr=p,nc=p)
  if (is.null(choix)) choix<-"eta"
  switch(choix,
         "eta"={
           if (forceEta){
             for (j in (1:p)) { for (k in (1:p)){
               if((H[j]+H[k])!=1 ) eta[j,k]<-eta[j,k]/(1-H[j]-H[k])
               
             }}
           }
         }, 
         "causal"={
           eta<-matrix(0,nc=p,nr=p)
           for (j in (1:p)) 
             for (k in (1:p)){
               Hjk<-H[j]+H[k]
               if (Hjk!=1){ eta[j,k]<- rho[j,k]*(cos(pi*H[j])-cos(pi*H[k]))/(cos(pi*H[j])+cos(pi*H[k]))	}
               else{ eta[j,k]<-rho[j,k]*2/pi/tan(pi*H[j])}
             }
           diag(eta)<-0
         },
         "wb"={
           eta<-matrix(0,nc=p,nr=p)
         })	
  
  ## computation of the matrix G (existence proposition)
  for (j in (1:p)){for (k in (1:p)){
    Hjk<-H[j]+H[k]
    if (Hjk!=1){
      pre<-rho[j,k]*sin(pi/2*Hjk) 
      pim<- -eta[j,k]*cos(pi/2*Hjk)
    }
    else{
      pre<-rho[j,k];pim<- -pi/2*eta[j,k]
    }
    G[j,k]<-complex(real=pre,imaginary=pim)
    G[j,k]<-gamma(H[j]+H[k]+1)*G[j,k]
  }}
  ltmp<-sum(eigen(G)$values<0)
  #print(H)
  #ltmp<-existMFBM(H=H,rho=rho,eta=eta,forceEta=forceEta,choix=choix)  
  txt<-paste("cannot be simulated for these parameters :",ltmp," negative eigenvalues")
  if (print & (ltmp==0)) cat('No negative eigenvalues for G \n')
  if (ltmp>0) stop(txt)
  
  xlogx<-function(h){
    res<-NULL
    for (e in h) { if (e==0) res<-c(res,0) else res<-c(res,h*log(abs(h)))}
    res
  }
  
  gammauv<-function(h,Hj,Hk,rhojk,etajk){
    #pour h>=1
    Hjk<-Hj+Hk
    wjk<-function(h,rhojk,etajk){
      tmp2<-NULL
      for (e in h){
        if (Hjk!=1) {
          tmp<-(rhojk-etajk*sign(e))*abs(e)^Hjk 
        }	else { tmp<-rhojk*abs(e)-etajk*xlogx(e)}
        tmp2<-c(tmp2,tmp)
      }
      tmp2
    }
    .5*(wjk(h-1,rhojk,etajk)-2*wjk(h,rhojk,etajk)+wjk(h+1,rhojk,etajk))
  }
  cuvj <- function(u,v, m){
    ## provides the first line of the matrix C_m(H_u,H_v) de Wood et Chan version p>1 !
    Hu<-H[u]
    Hv<-H[v]
    Huv<-Hu+Hv
    z0<-rho[u,v]
    z1<-gammauv(1:(m/2-1),Hu,Hv,rho[u,v],eta[u,v])
    z2<-gammauv(m/2,Hu,Hv,rho[u,v],eta[u,v])+gammauv(m/2,Hv,Hu,rho[v,u],eta[v,u])
    z2<-z2/2
    z3<-gammauv(m-(m/2+1):(m-1),Hv,Hu,rho[v,u],eta[v,u])	
    c(z0,z1,z2,z3)
  }
  m<-2^(trunc(log(n)/log(2))+2)
  tab<-B<-array(0,dim=c(p,p,m))
  for (u in 1:p){ for (v in 1:p){
    vpCuv <- cuvj( u, v,m)
    vpCuv <- (fft(c(vpCuv), inverse = FALSE))
    vpCuv<-(vpCuv)
    tab[u,v,]<-vpCuv
  }}
  res<-list()
  for (iSamp in (1:nbSamp)){
    cat('\r i=',iSamp)
    W<-matrix(0,nr=m,ncol=p)
    ## Simulations of U and V and storage in an array
    Z<-matrix(0,nr=m,nc=p)
    Z[1,]<-rnorm(p)/sqrt(m)
    Z[m/2+1,]<-rnorm(p)/sqrt(m)
    
    trace<-NULL
    for (j in (0:(m-1))){
      if ((j>0) & (j<=(m/2-1))){
        U<-rnorm(p); V<-rnorm(p)
        Z[j+1,]<-complex(real=U,imaginary=V) /sqrt(2*m)
        Z[m-j+1,]<-complex(real=U,imaginary=(-V)) /sqrt(2*m)
      } 
      A<-tab[,,j+1]
      tmp<-eigen(A)
      vpA<-Re(tmp$values)
      vecpA<-	tmp$vec
      vecpAinv<-solve((tmp$vec))
      vpAPlus<-apply(cbind(vpA,rep(0,length(vpA))),1,max)
      vpANeg<-apply(cbind(-vpA,rep(0,length(vpA))),1,max)
      trace<-c(trace,sum(vpANeg)/sum(vpA))
      B[,,j+1]<-vecpA %*% diag(sqrt((vpAPlus))) %*% vecpAinv
      W[j+1,]<-B[,,j+1] %*% Z[j+1,]
    }
    X<-Re(mvfft(W,inverse=F)[1:n,])
    vFBM<-apply(X,2,cumsum)
    if (print) cat('  Mean trace:',mean(trace),'\n')
    res[[iSamp]]<-vFBM
  }
  if (nbSamp==1) {  res2<- res[[1]]} else res2<-res	
  res2
  
}

mfBB <- function(n=50,end=1, H=c(.2,.2),sig=c(1,1),rhomat=matrix(c(1,.5,.5,1),nc=2),
                 etamat=matrix(c(0,.5,-.5,0),nc=2),nbSamp=1){ 
  # simulation of mfBB over interval (0,end] at discrete times end/n, 2*end/n, ..., end.
  # H: Hurst parameter vector of dimension p
  # sig: Vector of standard deviations (at time 1)
  # rhomat: matrix  with pairwise correlations (at time 1)
  # etamat: matrix with asymmetry parameters
  # nbSamp: No. of iterations. If =1, output is one (n,p)-matrix,
  #   if >1, output is a list of (n,p)-matrices.
  # Based on: https://sites.google.com/site/homepagejfc/software
  p<-length(H)
  # Check existence
  ltmp<-existMFBB(H=H,rhomat=rhomat,eta=etamat)  
  txt<-paste("Does not exist :",ltmp," negative eigenvalues")
  if (ltmp>0) stop(txt)
  xlogx<-function(h){
    res<-NULL
    for (e in h) { if (e==0) res<-c(res,0) else res<-c(res,h*log(abs(h)))}
    res
  }
  # autocovariance function of increments for vector h with lags >=1
  gammauv<-function(h,sj,sk,Hj,Hk,rhojk,etajk){
    Hjk<-Hj+Hk
    wjk<-function(h,rhojk,etajk){
      tmp2<-NULL
      for (e in h){
        if (Hjk!=1) {
          tmp<-(rhojk-etajk*sign(e))*abs(e)^Hjk 
        }	else { tmp<-rhojk*abs(e)+etajk*xlogx(e)}
        tmp2<-c(tmp2,tmp)
      }
      tmp2
    }
    sj*sk*.5*(wjk(h+1,rhojk,etajk)-2*wjk(h,rhojk,etajk)+wjk(h-1,rhojk,etajk))
  }
  # Simulation of fractional Gaussian noise at 0, ..., n-1
  cuvj <- function(u,v,m){
    su <- sig[u]
    sv <- sig[v]
    Hu<-H[u]
    Hv<-H[v]
    Huv<-Hu+Hv
    z0<-sig[u]*sig[v]*rhomat[u,v]#modified compared to Coeurjolly by adding factors sig[u]*sig[v]
    z1<-gammauv(1:(m/2-1),su,sv,Hu,Hv,rhomat[u,v],etamat[u,v])
    z2<-gammauv(m/2,su,sv,Hu,Hv,rhomat[u,v],etamat[u,v])+gammauv(m/2,su,sv,Hv,Hu,rhomat[v,u],etamat[v,u])
    z2<-z2/2
    z3<-gammauv(m-(m/2+1):(m-1),su,sv,Hv,Hu,rhomat[v,u],etamat[v,u])	
    c(z0,z1,z2,z3)
  }
  m<-2^(trunc(log(n)/log(2))+2)
  tab<-B<-array(0,dim=c(p,p,m))
  for (u in 1:p){ for (v in 1:p){
    vpCuv <- cuvj( u, v,m)
    vpCuv <- (fft(c(vpCuv), inverse = FALSE))
    vpCuv<-(vpCuv)
    tab[u,v,]<-vpCuv
  }}
  res<-list()
  for (iSamp in (1:nbSamp)){
    W<-matrix(0,nr=m,ncol=p)
    Z<-matrix(0,nr=m,nc=p)
    Z[1,]<-rnorm(p)/sqrt(m)
    Z[m/2+1,]<-rnorm(p)/sqrt(m)
    trace<-NULL
    for (j in (0:(m-1))){
      if ((j>0) & (j<=(m/2-1))){
        U<-rnorm(p); V<-rnorm(p)
        Z[j+1,]<-complex(real=U,imaginary=V) /sqrt(2*m)
        Z[m-j+1,]<-complex(real=U,imaginary=(-V)) /sqrt(2*m)
      } 
      A<-tab[,,j+1]
      tmp<-eigen(A)
      vpA<-Re(tmp$values)
      vecpA<-	tmp$vec
      vecpAinv<-solve((tmp$vec))
      vpAPlus<-apply(cbind(vpA,rep(0,length(vpA))),1,max)
      vpANeg<-apply(cbind(-vpA,rep(0,length(vpA))),1,max)
      trace<-c(trace,sum(vpANeg)/sum(vpA))
      B[,,j+1]<-vecpA %*% diag(sqrt((vpAPlus))) %*% vecpAinv
      W[j+1,]<-B[,,j+1] %*% Z[j+1,]
    }
    X<-Re(mvfft(W,inverse=F)[1:n,])
    # mfBB as cumulative sum of fractional Gaussian noise
    vFBM<-apply(X,2,cumsum)
    # Self-similarity to go from 1, ..., n to times end/n, ..., end.
    for (u in 1:p){
      vFBM[,u] <- rep(end/n,n)^H[u]*vFBM[,u]
    }
    res[[iSamp]]<-vFBM
  }
  if (nbSamp==1) {  res2<- res[[1]]} else res2<-res
  return(res2)
}

# Copied from: https://sites.google.com/site/homepagejfc/software
filt<-function(nm="i2"){
  ## gives the coefficients of an increments type filter, a Daublet, Symmlet or Coiflet filter
  fact<-function(n) ifelse(n==0,1,prod(1:n))
  cnk<-function(n,k) fact(n)/fact(k)/fact(n-k)
  
  if (!is.na(o<-as.numeric(strsplit(nm,"i")[[1]][2]))) {
    l<-o+1
    a<-rep(0,l)
    for (k in 0:o){
      a[k+1]<-cnk(o,k)*(-1)^(k)
    }
  }
  else {
    require(wmtsa)
    a<-wavDaubechies(nm)$wav
  }
  a
}
# Copied from: https://sites.google.com/site/homepagejfc/software
dilatation <- function(a=c(1,-2,1), m=2){
  ## provides the dilated version of a filter a
  if(m>1){
    la <- length(a)
    am <- rep(0, m * la - 1)
    am[seq(1, m * la - 1, by = m)] <- a
    am 
  } else a
}
# Copied from: https://sites.google.com/site/homepagejfc/software
dvFBM<-function(fbm,nma="i2",M1=1,M2=5,method=c("ST","Q","TM","B1-ST","B1-Q","B1-TM","B0-ST","B0-Q","B0-TM"),par=list(),llplot=FALSE){
  
  if (missing(fbm)) stop("Missing data")
  
  listDaub<-c("d2", "d4", "d6", "d8","d10", "d12", "d14", "d16", "d18", "d20","s2","s4", "s6", "s8", "s10","s12", "s14", "s16", "s18", "s20", "l2","l4", "l6", "l14", "l18", "l20","c6", "c12", "c18", "c24", "c30")
  l<-strsplit(nma,"")[[1]][1]
  if (l=="i") {
    if (length(strsplit(nma,"i")[[1]])>2) stop("Bad entry name for the filter ")
    if (is.na(as.numeric(strsplit(nma,"i")[[1]][2]))) stop("Bad entry name for the filter")
    a<-filt(nma)
  }else{
    if (nma %in% listDaub) {
      a<-filt(nma)
    }else stop("Bad entry name for the filter")
  }
  
  if (method %in% c("Q","B0-Q","B1-Q")){
    if (is.null(par$vecp) | is.null(par$vecc)) stop('"par=list(vecp=...,vecc=...) is needed for methods "Q","B0-Q","B1-Q"')
  }
  if (method %in% c("TM","B0-TM","B1-TM")){
    if (is.null(par$beta1) | is.null(par$beta2)) stop('par=list(beta1=...,beta2=...) is needed for methods "TM","B0-TM","B1-TM"')
  }
  
  if (!(method %in% c("ST","Q","TM","B1-ST","B1-Q","B1-TM","B0-ST","B0-Q","B0-TM"))) stop('Method should be one of "ST","Q","TM","B1-ST","B1-Q","B1-TM","B0-ST","B0-Q","B0-TM"')
  
  
  l<-length(a)-1
  n<-length(fbm)
  Unam<-NULL
  Unam<-switch(method,
               "ST"={
                 for (m in M1:M2){
                   am<-dilatation(a,m)
                   Vam <- filter(fbm, am, sides = 1)
                   Vam <- Vam[!is.na(Vam)]
                   Unam<-c(Unam,mean(Vam^2))
                 }
                 Unam
                 
               },
               "Q"={
                 vecp<-par$vecp
                 vecc<-par$vecc
                 for (m in M1:M2) {
                   am<-dilatation(a,m)
                   Vam <- filter(fbm,dilatation(a,m),sides=1)
                   Vam <- Vam[!is.na(Vam)]
                   Unam <- c(Unam, sum(vecc * quantile(Vam^2,vecp)))
                 }
                 Unam
               },
               "TM"={
                 beta1<-par$beta1
                 beta2<-par$beta2
                 for (m in M1:M2) {
                   am<-dilatation(a,m)
                   Vam <- filter(fbm, am, sides = 1)
                   Vam <- Vam[ - (1:(m*l))]
                   tmp<-sort(Vam^2)
                   nn<-length(Vam)
                   Unam<-c(Unam,mean(tmp[(trunc(nn*beta1)+1):(nn-trunc(nn*beta2))]))
                 }
                 Unam
               },
               "B1-ST"={
                 for (m in M1:M2) {
                   Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
                   Va2m <- Va2m[!is.na(Va2m)]
                   Vam <- filter(fbm,dilatation(a,m),sides=1)
                   Vam <- Vam[!is.na(Vam)]
                   Unam <- c(Unam,abs(mean(Va2m^2)-mean(Vam^2)))
                 }
                 Unam
               },
               "B1-Q"={
                 vecp<-par$vecp
                 vecc<-par$vecc
                 for (m in M1:M2) {
                   Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
                   Va2m <- Va2m[!is.na(Va2m)]
                   Vam <- filter(fbm,dilatation(a,m),sides=1)
                   Vam <- Vam[!is.na(Vam)]
                   Unam <- c(Unam, abs(sum(vecc * quantile(Va2m^2,vecp)) - sum(vecc * quantile(Vam^2,vecp))))
                 }
                 Unam
               },
               "B1-TM"={
                 beta1<-par$beta1
                 beta2<-par$beta2
                 for (m in M1:M2) {
                   Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
                   Va2m <- Va2m[!is.na(Va2m)]
                   Vam <- filter(fbm,dilatation(a,m),sides=1)
                   Vam <- Vam[!is.na(Vam)]
                   tmp<-sort(Vam^2);nn<-length(Vam)
                   tmp2<-sort(Va2m^2);nn2<-length(Va2m)
                   Unam<-c(Unam, abs(mean(tmp2[(trunc(nn2*beta1)+1):(nn2-trunc(nn2*beta2))]) -mean(tmp[(trunc(nn*beta1)+1):(nn-trunc(nn*beta2))])))
                 }
                 Unam
               },
               "B0-ST"={
                 for (m in M1:M2) {
                   Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
                   Va2m <- Va2m[!is.na(Va2m)]
                   Vam <- filter(fbm,dilatation(a,m),sides=1)
                   Vam <- Vam[!is.na(Vam)]
                   Unam <- c(Unam,abs(mean(Va2m^2)/(2*m)-mean(Vam^2)/m))
                 }
                 Unam
               },
               "B0-Q"={
                 vecp<-par$vecp
                 vecc<-par$vecc
                 for (m in M1:M2) {
                   Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
                   Va2m <- Va2m[!is.na(Va2m)]
                   Vam <- filter(fbm,dilatation(a,m),sides=1)
                   Vam <- Vam[!is.na(Vam)]
                   Unam <- c(Unam, abs(sum(vecc * quantile(Va2m^2,vecp))/(2*m) - sum(vecc * quantile(Vam^2,vecp))/m))
                 }
                 Unam
               },
               "B0-TM"={
                 beta1<-par$beta1
                 beta2<-par$beta2
                 for (m in M1:M2) {
                   Va2m <- filter(fbm,dilatation(a,2*m),sides=1)
                   Va2m <- Va2m[!is.na(Va2m)]
                   Vam <- filter(fbm,dilatation(a,m),sides=1)
                   Vam <- Vam[!is.na(Vam)]
                   tmp<-sort(Vam^2);nn<-length(Vam)
                   tmp2<-sort(Va2m^2);nn2<-length(Va2m)
                   Unam<-c(Unam, abs(mean(tmp2[(trunc(nn2*beta1)+1):(nn2-trunc(nn2*beta2))])/(2*m) -mean(tmp[(trunc(nn*beta1)+1):(nn-trunc(nn*beta2))])/m ))
                 }
                 Unam
               }
  )
  reg<-lm(log(Unam)~log(M1:M2))
  opt<-rev(reg$coeff)
  
  if (method %in% c("B0-ST","B0-Q","B0-TM")) Hest<-(opt[1]+1)/2
  else Hest<-opt[1]/2
  
  if (llplot) { plot(log(Unam)~log(M1:M2));abline(reg)}
  Hest
}

pia1a2H<-function(a1,a2,H,j){
  l1<-length(a1)-1
  l2<-length(a2)-1
  s<-0
  for (q in (0:l1)){for (r in (0:l2)){
    s<-s+a1[q+1]*a2[r+1]*(abs(q-r+j))^(2*H)
  }}
  -.5*s
}
covEmp<-function(xi,xj,h){
  ## computes the empirical covariance of xi and xj at lag h
  n<-length(xi)
  mean(xi[1:(n-h)]* xj[(h+1):n])
}
# Copied from: https://sites.google.com/site/homepagejfc/software
estMFBM<-function(x,nma="i2",M1=1,M2=5,w=c(1,0,0),forceH=FALSE,naiveRhoEta=FALSE){
  ## Parameters estimation of a MFBM using discrete variations techniques.
  ## x= a sample path of a MFBM.
  ## w= weight vector = (wv,wc,wd) wc=1 use of instanteneous covariances, wd=1 use of g_ij(+-m)
  ## Notation follow the ones in the related reference.
  ##
  ## Reference. Amblard, Coeurjolly (2011)
  ## Identification of the multivariate fractional Brownian motion. 
  ## IEEE Transactions on Signal Processing 59(11), 5152–5168.
  ##
  ##
  p<-ncol(x)
  n<-nrow(x)
  wv<-w[1];wc<-w[2];wd<-w[3]
  lambda<-4*wv+(wc+wd)*(p-2)
  lambdap<-lambda+(wc+wd)*p
  J<-matrix(1,nr=p,nc=p)
  Ainv<-1/p*(1/lambdap-1/lambda)*J
  diag(Ainv)<-diag(Ainv)+1/lambda
  a<-filt(nma)
  ell<-length(a)-1
  b<-bv<-NULL 
  
  #if (M1>=M2) stop("M1 should be < M2")
  #if ((wd==1) & (M1<= (length(a)-1))) stop("If wd!=0 M1 should be > ell=length(a)-1")
  LM<-log(M1:M2);LM<-LM-mean(LM)
  for (i in (1:p)){
    vi<-scij<-sdij<-NULL
    for (m in M1:M2) {
      am <- dilatation(a, m)
      xiam <- filter(x[,i], am, sides = 1)
      xiam <- xiam[!is.na(xiam)]	
      vi<-c(vi,log(mean(xiam^2)))
      cij<-dij<-0
      for (j in ((1:p)[-i]) ){
        xjam<-filter(x[,j], am, sides = 1)
        xjam <- xjam[!is.na(xjam)]
        nn<-length(xiam)
        cij<-cij+ log(abs(mean(xiam*xjam)))
        #if(M1>=(length(a)-1)){
        cijm<-covEmp(xiam,xjam,m*ell)
        cjim<-covEmp(xjam,xiam,m*ell)
        #tmp1<-mean(xiam[1:(nn-m)]*xjam[(m+1):nn])
        #tmp2<-mean(xjam[1:(nn-m)]*xiam[(m+1):nn])
        dij<-dij+ log(abs(cijm-cjim))
        #} else dij<-0
      }
      scij<-c(scij,cij)
      sdij<-c(sdij,dij)
    }
    b<-c(b,t(LM)%*%(2*wv*vi+(wc*scij+wd*sdij) ) / (t(LM)%*%LM))
    bv<-c(bv,t(LM)%*%(vi) / (t(LM)%*%LM))
  }
  Hest<-Ainv %*% t(t(b))
  Hestv<-bv/2
  
  sEst<-NULL
  
  if (naiveRhoEta==TRUE) { M1<-M2<-1}
  
  if(forceH==TRUE) {Hest2<-as.vector(Hestv)}
  else Hest2<-as.vector(Hest)
  
  eps<-.001;Hest2[Hest2<0]<-eps;Hest2[Hest2>1]<-1-eps
  for (i in (1:p)){
    vi<-NULL
    for (m in M1:M2) {
      am <- dilatation(a, m)
      xiam <- filter(x[,i], am, sides = 1)
      xiam <- xiam[!is.na(xiam)]	
      vi<-c(vi,log(mean(xiam^2)))
      
    }
    
    sEst[i]<- exp(mean(vi-2*Hest2[i]*log(M1:M2)))/pia1a2H(a,a,Hest2[i],0)
  }
  
  
  rhoEst<-etaEst<-matrix(0,nr=p,nc=p)
  
  #if (wc==1){
  for (i in (1:(p-1))){ for (j in (i+1):p) {
    cij<-NULL
    for (m in M1:M2) {
      am <- dilatation(a, m)
      xiam <- filter(x[,i], am, sides = 1)
      xiam <- xiam[!is.na(xiam)]	
      xjam <- filter(x[,j], am, sides = 1)
      xjam <- xjam[!is.na(xjam)]
      cij<-c(cij,log(abs(  mean(xiam*xjam)) ) )
      
    }
    muijest<-mean(cij-(Hest2[i]+Hest2[j])*log(M1:M2))
    #print(muijest)
    #print(c(pia1a2H(a,a,(Hest[i]+Hest[j])/2,0),pia1a2H(a,a,(Hest[i]+Hest[i])/2,0),pia1a2H(a,a,(Hest[j]+Hest[j])/2,0),sqrt(sEst[i]*sEst[j])))
    if (M1==M2){
      rhoEst[i,j]<- mean(xiam*xjam)/sd(xiam)/sd(xjam)*sqrt(pia1a2H(a,a,Hest2[i],0)*pia1a2H(a,a,Hest2[j],0) )/pia1a2H(a,a,(Hest2[i]+Hest2[j])/2,0)
    } else { rhoEst[i,j]<-exp(muijest)/pia1a2H(a,a,(Hest2[i]+Hest2[j])/2,0)/sqrt(sEst[i]*sEst[j])}
    
    #print(rhoEst)
  }}
  #}
  #print(Hest2)
  #if (wd==1){
  for (i in (1:(p-1))){ for (j in (i+1):p) {
    dij<-NULL
    for (m in M1:M2) {
      am <- dilatation(a, m)
      xiam <- filter(x[,i], am, sides = 1)
      xiam <- xiam[!is.na(xiam)]	
      xjam <- filter(x[,j], am, sides = 1)
      xjam <- xjam[!is.na(xjam)]
      cijm<-covEmp(xiam,xjam,m*ell)
      cjim<-covEmp(xjam,xiam,m*ell)
      #cat("i,j",i,j,"\n");print(c(cijm,cjim))
      #print(c(mean(xiam^2),mean(xjam^2)))
      if (M1==M2) dij<-c(dij,abs(cijm-cjim))
      else dij<-c(dij, log(abs(cijm-cjim)))
      
    }
    nuijest<-mean(dij-(Hest2[i]+Hest2[j])*log(M1:M2))
    #print(c(pia1a2H(a,a,(Hest[i]+Hest[j])/2,0),pia1a2H(a,a,(Hest[i]+Hest[i])/2,0),pia1a2H(a,a,(Hest[j]+Hest[j])/2,0),sqrt(sEst[i]*sEst[j])))
    #etaEst[i,j]<-exp(nuijest)/denom4eta(a,(Hest[i]+Hest[j])/2)/sqrt(sEst[i]*sEst[j])
    if (M1==M2) {
      etaEst[i,j]<-dij/sd(xiam)/sd(xjam) *sqrt(pia1a2H(a,a,Hest2[i],0)*pia1a2H(a,a,Hest2[j],0) )  /abs(pia1a2H(a,a,(Hest2[i]+Hest2[j])/2,ell))/2 
      #print(dij);print(abs(pia1a2H(a,a,(Hest2[i]+Hest2[j])/2,ell)));print(sqrt(sEst[i]*sEst[j]))
      #print(etaEst[i,j])
      #print(dij/2/sqrt(mean(xiam^2)*mean(xjam^2))->poi1 );print(poi2<-pia1a2H(a,a,(Hest2[i]+Hest2[j])/2,ell)/sqrt(pia1a2H(a,a,Hest2[i],0)*pia1a2H(a,a,Hest2[j],0)))
      #print(poi1/poi2)
    }
    else etaEst[i,j]<-exp(nuijest)/abs(pia1a2H(a,a,(Hest2[i]+Hest2[j])/2,ell))/2 /sqrt(sEst[i]*sEst[j])
    #print(rhoEst)
  }}
  # }
  Hest[Hest<0]<-eps;Hest[Hest>1]<-1-eps
  
  list(H=Hest,s=sEst,rho=rhoEst,eta=etaEst)
}