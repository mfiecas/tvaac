library(nlme)
library(MASS)

tpower1<-function(x,t,p=1){(x-t)^p*(x>t)}

#generates truncated polynomial basis.
generateTruncatedBasis1<-function(n, knots,p=1){
  C<-Z<-NULL
  n.knots<-length(knots)
  vec<-1:n
  for (i in 1:(n.knots) ){
    Z<-cbind(Z,tpower1(vec, knots[i],p=p))
  }
  Z<-Z/max(Z)
  if (p>=1){
    for (j in 1:p){
      C<-cbind(C, ((1:n)/n)^j )
    }
  }
  out<-cbind(C,Z)
  return(out)
}


#computes time-varying activation and connectivity
tvaac<-function(y,X,W=NULL,AR=1,tv=NULL, off.diag=NULL, n.knots=NULL,t.knots=NULL,power=1){
  #y=(n x 2) matrix
  #X=(n x p) matrix corresponding to the time-varying effects, currenty supports p=1 or 2
  #W=(n x q) matrix of adjusted covariates 
  #tv="both" if time-varying mean and variance
  #  ="mean" if time-varying mean only
  #  ="var"  if time-varying correlation only
  #  ="none" if no temporal dynamics are applied.
  #AR=AR parameters, cuurently supports 1,2 or 3. >3 will be suppressed to 3.
  #n.knots= number of knots
  #t.knots= pre-specified knots
  #power=degree of truncated polynomial basis
  
  X<-as.matrix(X)
  n<-nrow(X)
  p<-ncol(X)
  
  if (is.null(off.diag)){ 
    off.diag<-"constant"; 
    print("off.diag set to constant")}
  if (is.null(tv) || !(tv%in%c("both","mean", "var","none")) ){
    tv<-"both"
    print("tv argument set to both")
  }
  
  if (AR>3){ 
    print(paste0("AR=",AR," is not supported and reduced to AR(3)" ));
    AR=3
  }
  
  
  if (is.null(W)){ W<-matrix(rep(1,n)) }
  else{
    q<-ncol(W)
    bool=0
    for (j in 1:q){if (isTRUE(all.equal(W[,j],rep(1, n)) )) {bool<-1}}
    if (bool==0){ W<-cbind(W, rep(1,n)) }
  }
  y1<-y[,1]
  y2<-y[,2]
  
  if (tv %in% c("both", "mean", "var")){
    if (!is.null(t.knots)){
      if (!is.null(n.knots)){ print("n.knots ignored as t.knots is specified") }
    }
    else{
      if (is.null(n.knots)){
        print("number of knots set to 10")
        n.knots=10
      }
      t.knots<-seq(0, n, by=n/n.knots)[-1]
      t.knots<-t.knots[-length(t.knots)]
    }
    
    Z<-cbind(generateTruncatedBasis1(n, t.knots,p=power))
  }
  
  betahat.mat<-betahat.mat2<-matrix(NA, n, p)
  
  if (tv%in%c("both", "mean")){
    if (p==1){
      all<-rep(1,n) 
      Zstim1<-Z
      for (j in 1:ncol(Z)){ Zstim1[,j]<-Z[,j]*X }
      mod1<-lme(y1~-1+X+W, random=list(all=pdIdent(~Zstim1-1)))
      mod2<-lme(y2~-1+X+W, random=list(all=pdIdent(~Zstim1-1))) 
      
      betahat.mat[,1]=mod1$coefficients$fixed[1]+Z%*%c(mod1$coefficients$random[[1]])
      betahat.mat2[,1]=mod2$coefficients$fixed[1]+Z%*%c(mod2$coefficients$random[[1]])
      
    }  else if (p==2){
      all<-all2<-rep(1,n)
      X1<-X[,1]
      X2<-X[,2]
      Zstim1<-Zstim2<-Z
      Zstim0a<-Zstim0b<-Z
      for (j in 1:ncol(Z)){
        Zstim0a[,j]<-Z[,j]*(X1+X2)
        Zstim0b[,j]<-Z[,j]*(X1-X2)
        Zstim1[,j]<-Z[,j]*X1
        Zstim2[,j]<-Z[,j]*X2
      }
      
      mod1a<-lme(y1~-1+X+W, random=list(all=pdIdent(~Zstim0a-1), all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1)))
      mod1b<-lme(y1~-1+X+W, random=list(all=pdIdent(~Zstim0b-1), all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1)))
      
      if (logLik(mod1a)>logLik(mod1b)){ 
        mod1<-mod1a; rm(mod1b); rm(mod1a)
        a<- 1
      }
      else{ 
        mod1<-mod1b; rm(mod1a); rm(mod1b)
        a<- -1
      }
      
      mod2a<-lme(y2~-1+X+W, random=list(all=pdIdent(~Zstim0a-1), all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1)))
      mod2b<-lme(y2~-1+X+W, random=list(all=pdIdent(~Zstim0b-1), all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1)))
      if (logLik(mod2a)>logLik(mod2b)){ 
        mod2<-mod2a; rm(mod2a); rm(mod2b)
        b<- 1
      }
      else{ 
        mod2<-mod2b; rm(mod2a); rm(mod2b) 
        b<- -1
      }
      betahat.mat[,1]<- mod1$coefficients$fixed[1]+Z%*%c(mod1$coefficients$random[[1]]+mod1$coefficients$random[[2]])
      betahat.mat[,2]<- mod1$coefficients$fixed[2]+Z%*%c(a*mod1$coefficients$random[[1]]+mod1$coefficients$random[[3]])
      betahat.mat2[,1]<- mod2$coefficients$fixed[1]+Z%*%c(mod2$coefficients$random[[1]]+mod2$coefficients$random[[2]])
      betahat.mat2[,2]<- mod2$coefficients$fixed[2]+Z%*%c(b*mod2$coefficients$random[[1]]+mod2$coefficients$random[[3]])
    } 
  }
  else{
    mod1<-lm(y1~-1+X+W)
    mod2<-lm(y2~-1+X+W)
    for (k in 1:p){
      betahat.mat[,k]<- mod1$coefficients[k]
      betahat.mat2[,k]<-mod2$coefficients[k]
    }
  }
  
  resid1<-y1-predict(mod1)
  resid2<-y2-predict(mod2)
  
  ind1<-AR+1
  ind2<-n
  
  resid1.y<-resid1[ind1:ind2]
  resid2.y<-resid2[ind1:ind2]
  
  resid1.X<-NULL
  resid2.X<-NULL
  
  for (ar in 1:AR){
    index<-(ind1-ar):(ind2-ar)
    resid1.X<-cbind(resid1.X, resid1[index])
    resid2.X<-cbind(resid2.X, resid2[index])
    
    assign(paste0("all",ar), rep(1,n-AR))
  }
  
  
  if (tv%in%c("both", "var")){
    all<-rep(1, n-AR)
    Z1<-generateTruncatedBasis1(n-AR, AR+t.knots, p=power)
    if (AR==1){
      index<-(ind1-1):(ind2-1)
      Zstim1<-ZZstim1<-matrix(NA, n-AR, ncol(Z1))
      for (j in 1:ncol(Z1)){ 
        Zstim1[,j]<-resid1[index]*Z1[,j]   
        ZZstim1[,j]<-resid2[index]*Z1[,j]
      }
      if (off.diag=="constant"){
        mod1.AR<-lme(resid1.y~-1+resid1.X+resid2.X, random=list(all=pdIdent(~Zstim1-1)))
        mod2.AR<-lme(resid2.y~-1+resid2.X+resid1.X, random=list(all=pdIdent(~ZZstim1-1)))
      }
      else if (off.diag=="tv"){
        mod1.AR<-lme(resid1.y~-1+resid1.X+resid2.X, random=list(all=pdIdent(~Zstim1-1), all=pdIdent(~ZZstim1-1)))
        mod2.AR<-lme(resid2.y~-1+resid2.X+resid1.X, random=list(all=pdIdent(~ZZstim1-1),all=pdIdent(~Zstim1-1)))
      }
      
    }
    else if (AR==2){
      index<-(ind1-1):(ind2-1)
      index2<-(ind1-2):(ind2-2)
      
      Zstim1<-ZZstim1<-Zstim2<-ZZstim2<-matrix(NA, n-AR, ncol(Z1))
      for (j in 1:ncol(Z1)){ 
        Zstim1[,j]<-resid1[index]*Z1[,j]   
        ZZstim1[,j]<-resid2[index]*Z1[,j]
        Zstim2[,j]<-resid1[index2]*Z1[,j]   
        ZZstim2[,j]<-resid2[index2]*Z1[,j]
      }
      if (off.diag=="constant"){
        mod1.AR<-lme(resid1.y~-1+resid1.X+resid2.X, random=list(all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1) ))
        mod2.AR<-lme(resid2.y~-1+resid2.X+resid1.X, random=list(all=pdIdent(~ZZstim1-1), all=pdIdent(~ZZstim2-1) ))
      }
      else if(off.diag=="tv"){
        mod1.AR<-lme(resid1.y~-1+resid1.X+resid2.X, random=list(all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1), all=pdIdent(~ZZstim1-1),all=pdIdent(~ZZstim2-1)  ))
        mod2.AR<-lme(resid2.y~-1+resid2.X+resid1.X, random=list(all=pdIdent(~ZZstim1-1), all=pdIdent(~ZZstim2-1),all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1) ))
      }
      
    }
    else {
      index<-(ind1-1):(ind2-1)
      index2<-(ind1-2):(ind2-2)
      index3<-(ind1-3):(ind2-3)
      
      Zstim1<-ZZstim1<-Zstim2<-ZZstim2<-Zstim3<-ZZstim3<-matrix(NA, n-AR, ncol(Z1))
      for (j in 1:ncol(Z1)){ 
        Zstim1[,j]<-resid1[index]*Z1[,j]   
        ZZstim1[,j]<-resid2[index]*Z1[,j]
        Zstim2[,j]<-resid1[index2]*Z1[,j]   
        ZZstim2[,j]<-resid2[index2]*Z1[,j]
        Zstim3[,j]<-resid1[index3]*Z1[,j]   
        ZZstim3[,j]<-resid2[index3]*Z1[,j]
        
      }
      if (off.diag=="constant"){
        mod1.AR<-lme(resid1.y~-1+resid1.X+resid2.X, random=list(all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1), all=pdIdent(~Zstim3-1) ))
        mod2.AR<-lme(resid2.y~-1+resid2.X+resid1.X, random=list(all=pdIdent(~ZZstim1-1), all=pdIdent(~ZZstim2-1), all=pdIdent(~ZZstim3-1) ))
      }
      else if (off.diag=="tv"){
        mod1.AR<-lme(resid1.y~-1+resid1.X+resid2.X, random=list(all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1), all=pdIdent(~Zstim3-1),
                                                                all=pdIdent(~ZZstim1-1), all=pdIdent(~ZZstim2-1), all=pdIdent(~ZZstim3-1)))
        mod2.AR<-lme(resid2.y~-1+resid2.X+resid1.X, random=list(all=pdIdent(~ZZstim1-1), all=pdIdent(~ZZstim2-1), all=pdIdent(~ZZstim3-1),
                                                                all=pdIdent(~Zstim1-1), all=pdIdent(~Zstim2-1), all=pdIdent(~Zstim3-1)))
      }
    }
  }
  else{
    mod1.AR<-lm(resid1.y~-1+resid1.X+resid2.X)
    mod2.AR<-lm(resid2.y~-1+resid2.X+resid1.X)
  }
  
  
  
  delta1<-resid1.y-predict(mod1.AR)
  delta2<-resid2.y-predict(mod2.AR)
  
  Sigma.est<-cov(cbind(delta1,delta2))
  
  
  phihat1.mat<-phihat2.mat<-matrix(NA, n-AR, AR)
  phihat12.mat<-phihat21.mat<-matrix(NA, n-AR, AR)
  for (ar in 1:AR){
    if (tv%in%c("both", "var")){
      phihat1.mat[,ar]<-mod1.AR$coefficients$fixed[ar]+Z1%*%c(mod1.AR$coefficients$random[[ar]])
      phihat2.mat[,ar]<-mod2.AR$coefficients$fixed[ar]+Z1%*%c(mod2.AR$coefficients$random[[ar]])
      phihat12.mat[,ar]<-mod1.AR$coefficients$fixed[AR+ar]
      phihat21.mat[,ar]<-mod2.AR$coefficients$fixed[AR+ar]
      if (off.diag=="tv"){
        phihat12.mat[,ar]<-phihat12.mat[,ar]+Z1%*%c(mod1.AR$coefficients$random[[AR+ar]])
        phihat21.mat[,ar]<-phihat21.mat[,ar]+Z1%*%c(mod2.AR$coefficients$random[[AR+ar]])
      }
    }
    else{
      phihat1.mat[,ar]<-mod1.AR$coefficients[ar]
      phihat2.mat[,ar]<-mod2.AR$coefficients[ar]
      phihat12.mat[,ar]<-mod1.AR$coefficients[AR+ar]
      phihat21.mat[,ar]<-mod2.AR$coefficients[AR+ar]
    }
  }
  
  Blist<-list()
  for (ar in 1:AR){
    B<-array(0, c(2,2,n-AR))
    for (t in 1:(n-AR)){
      B[1,1,t]<-phihat1.mat[t,ar]
      B[1,2,t]<-phihat12.mat[t,ar]
      B[2,1,t]<-phihat21.mat[t,ar]
      B[2,2,t]<-phihat2.mat[t,ar]
    }
    Blist[[paste0("B",ar)]]<-B
  }
  
  tvcor.est = rep(NA,n-AR)
  tvcov.est = rep(NA,n-AR)
  for(t in 1:(n-AR)){
    Psi = matrix(0,2*AR,2*AR)
    ind<-1:2
    for (ar in 1:AR){
      Psi[1:2, ind]<-Blist[[ar]][,,t]
      if (AR>1 & ind[2]<nrow(Psi)){
        Psi[ind+2,ind]<-diag(2)
        ind<-ind+2
      }
    }
    r<-solve(diag(rep(1,(2*AR)^2))-kronecker(Psi,Psi))%*%c(kronecker(diag(AR),Sigma.est))
    r.mat<-matrix(r,2*AR,2*AR)
    tvcor.est[t] = cov2cor(r.mat)[1,2]
    tvcov.est[t] = r.mat[1,2]
  }
  
  return(
    list(
      y=y,
      X=X,
      W=W,
      power=power,
      off.diag=off.diag,
      resid=cbind(resid1, resid2),
      delta=cbind(delta1, delta2),
      betahat.mat=betahat.mat,
      betahat.mat2=betahat.mat2,
      phihat.mat1=phihat1.mat,
      phihat.mat2=phihat2.mat,
      phihat.mat12=phihat12.mat,
      phihat.mat21=phihat21.mat,
      tvcor.est=tvcor.est,
      tvcov.est=tvcov.est,
      B=Blist,
      AR=AR,
      mod1=mod1,
      mod2=mod2,
      mod1.AR=mod1.AR,
      mod2.AR=mod2.AR,
      t.knots=t.knots,
      tv=tv
    )
  )
}


#Bootstrapping tvaac for inference
boot.tvaac<-function(fit, sim.size=100, boot="parametric",seed=1131, pbar=TRUE){
  #fit: fit from tvmean.mv()
  #sim.size : # of bootstrap sample
  #method : "parametric" assumes normality of delta 
  #         "residual" uses resampling
  set.seed(seed)
  AR<-fit$AR
  Blist<-fit$B
  X<-fit$X
  W<-fit$W
  tv<-fit$tv
  off.diag<-fit$off.diag
  power<-fit$power
  t.knots<-fit$t.knots
  betahat.mat<-fit$betahat.mat
  betahat.mat2<-fit$betahat.mat2
  tvcor.est<-fit$tvcor.est
  p<-ncol(X)
  
  resid<-fit$resid
  mod1<-fit$mod1
  mod2<-fit$mod2
  pred.mod1<-predict(mod1)
  pred.mod2<-predict(mod2)
  delta<-fit$delta
  delta.mean<-apply(delta,2,mean)
  for (j in 1:nrow(delta)){
    delta[j,]<-delta[j,]-delta.mean
  }
  covDelta<-cov(delta)
  
  betahat.mat.boot<-array(NA, c(dim(betahat.mat),sim.size))
  betahat.mat2.boot<-array(NA, c(dim(betahat.mat2),sim.size))
  tvcor.boot<-matrix(NA, length(fit$tvcor.est), sim.size)
  tvcov.boot<-matrix(NA, length(fit$tvcov.est), sim.size)
  phihat.mat.boot<-array(NA, c(dim(fit$phihat.mat1), sim.size))
  phihat.mat2.boot<-array(NA, c(dim(fit$phihat.mat2), sim.size))
  
  if (pbar==TRUE) pb <- txtProgressBar(min = 0, max=sim.size, initial=0, char="-", style = 3)
  
  for (sim in 1:sim.size){
    tryCatch({
      epsilon<-resid
      for (j in (AR+1):nrow(resid)){
        est<-c(0,0)
        for (ar in 1:AR){
          B<-Blist[[ar]]
          est<-est+B[,,j-AR]%*%epsilon[j-ar,]
        }
        if (boot=="residual"){
          epsilon[j,]<-est+delta[sample(1:nrow(delta))[1],]
        }
        else if (boot=="parametric"){
          epsilon[j,]<-est+mvrnorm(1, c(0,0), covDelta)
        }
      }
      
      epsilon.mean<-apply(epsilon,2,mean)
      for (j in 1:nrow(epsilon)){
        epsilon[j,]<-epsilon[j,]-epsilon.mean
      }
      
      y=cbind(pred.mod1, pred.mod2)+epsilon
      fit.boot<-tvaac(y,X,W,tv = tv, AR=AR, n.knots = NULL, t.knots=t.knots, power = power, off.diag=off.diag)
      betahat.mat.boot[,,sim]<-fit.boot$betahat.mat
      betahat.mat2.boot[,,sim]<-fit.boot$betahat.mat2
      tvcor.boot[,sim]<-fit.boot$tvcor.est
      tvcov.boot[,sim]<-fit.boot$tvcov.est
      phihat.mat.boot[,,sim]<-fit.boot$phihat.mat1
      phihat.mat2.boot[,,sim]<-fit.boot$phihat.mat2
    }, 
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (pbar==TRUE){  setTxtProgressBar(pb, sim)  }
  }
  
  return(list(
    boot=boot,
    fit=fit,
    betahat.mat.boot=betahat.mat.boot,
    betahat.mat2.boot=betahat.mat2.boot,
    tvcor.boot=tvcor.boot,
    tvcov.boot=tvcov.boot,
    phihat.mat.boot=phihat.mat.boot,
    phihat.mat2.boot=phihat.mat2.boot
  ))
  
}


#Computes the confidence interval of interest
curve.conf.int<-function(fit.mat, method=NULL,conf.level=0.95){
  if (is.null(method) || !method%in%c("interval", "band")){ 
    print("method is set to confidence bands") 
    method="band" 
  }
  
  alpha<-1-conf.level
  ind.na<-which(!is.na(apply(fit.mat,2,sum)))
  fit.mat1<-fit.mat[,ind.na]
  sim.size<-ncol(fit.mat1)
  if (method=="interval"){
    conf.int<-t(apply(fit.mat1, 1, quantile, probs=c(alpha/2, 1-alpha/2)))
  }
  else if (method=="band"){
    dev<-rep(NA, sim.size)
    fit<-apply(fit.mat1,1,mean)
    for (j in 1:sim.size){
      dev[j]<-sum( (fit-fit.mat1[,j])^2)
    }
    ind<-which(dev<quantile(dev, 0.95))
    conf.int<-cbind(apply(fit.mat1[,ind],1,min), apply(fit.mat1[,ind],1,max))
  }
  return(conf.int)
}

