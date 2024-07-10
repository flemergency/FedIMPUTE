library(tidyverse)
library(pROC)
library(MASS)
library(glmnet)
library(stats)
library(Metrics)

lasso.fun=function(x,y, family, offset=NULL, weights=NULL, lambda=NULL,alpha=1){
  x = as.matrix(x)
  if(is.null(weights)){weights=rep(1,dim(x)[1])}
  if(is.null(lambda)){
    cvfit = glmnet::cv.glmnet(x, y, family=family, offset=offset, weights=weights, alpha=alpha)
    pred = predict(cvfit, newx=x, type = "response")
    coef.1se = coef(cvfit, s = "lambda.min")
  }else{
    fit = glmnet::glmnet(x, y, family=family, offset=offset, weights=weights,lambda=lambda, alpha=alpha)
    pred = predict(fit, newx = x, type = "response")
    coef.1se = coef(fit)
    CI = confint(fit, level = 0.95, parm = coef.1se)
  }
  best_cutoff = NA
  return(list(coef = as.numeric(coef.1se), cutoff = best_cutoff))
}


# translasso for calibration:
# default setting now: use target data to decide lambda
# modified: choose lambda using imputation model itself instead of downstream performances
# x.t and w.fit do not include intercept
Trans.fun=function(w.fit, x.t, y.t, n.s, p, lam.const=NULL, family, alpha=1, penalty=T, lam.TL.target=T){
  pred = x.t%*%w.fit
  mse.fed = mse(y.t,pred)
  if(penalty){
    lam.vec = 10^seq(-100,3,0.1)
    # choose the best w.fit using all target data
    w.fit.best = delta.fit.best = w.fit
    lam.best = 0
    MSE.best = Inf
    for(l in 1:length(lam.vec)){
      w.fit.temp <- w.fit*(abs(w.fit)>=lam.vec[l]*sqrt(2*log(p)/n.s))
      myoffset=x.t%*%w.fit.temp
      delta.fit <- as.numeric(glmnet::glmnet(x=x.t,y=y.t, offset=myoffset, lambda=lam.best*sqrt(2*log(p)/length(y.t)), family=family, alpha=alpha)$beta)
      delta.fit<-delta.fit*(abs(delta.fit)>=lam.vec[l]*sqrt(2*log(p)/length(y.t)))
      pred = x.t%*%(w.fit.temp+delta.fit)
  
      MSE.temp = mse(y.t,pred)
      
      if(MSE.temp < MSE.best){
        MSE.best = MSE.temp
        lam.best = lam.vec[l]
        w.fit.best = w.fit.temp
        delta.fit.best = delta.fit}
    }
  }
  else{
    myoffset=x.t%*%w.fit
    delta.fit <- as.numeric((glm(y.t~x.t, family = family, offset = myoffset))$coef)
    delta.fit.best = delta.fit[-1]
  }
  beta.fit <- w.fit + delta.fit.best
  beta.fit
}


### DAC main functions
# binomial:
iteration.fun2=function (dat.list, bini, kk.list) {
  K = length(dat.list)
  Uini.list = sapply(kk.list, function(kk) {
    U.fun(bini, dat = dat.list[[kk]])
  })
  Aini.list = lapply(1:K, function(kk) {
    A.fun(bini, dat = dat.list[[kk]])
  })
  Ahat.ini = Reduce("+", Aini.list)/K
  bhat.list = -ginv(Ahat.ini) %*% Uini.list + bini
  list(b.k = bhat.list, Ahat = Ahat.ini)
}
# gaussian:
iteration.fun3=function (dat.list, bini, kk.list) {
  K = length(dat.list)
  Uini.list = sapply(kk.list, function(kk) {
    U.fun.2(bini, dat = dat.list[[kk]])
  })
  Aini.list = lapply(1:K, function(kk) {
    A.fun.2(bini, dat = dat.list[[kk]])
  })
  Ahat.ini = Reduce("+", Aini.list)/K
  bhat.list = -ginv(Ahat.ini) %*% Uini.list + bini
  list(b.k = bhat.list, Ahat = Ahat.ini)
}
DCOS.FUN=function(dat.list, niter, ridge=T, lambda=NULL, family){
  K=length(dat.list)
  lambda.grid=10^seq(-100,3,0.1)
  N=sum(unlist(lapply(dat.list, function(xx) dim(xx)[1])))
  ### Round 1
  dat.list.00=dat.list[[1]]
  y1=dat.list.00[,1]
  x1=dat.list.00[,-1]
  # return a lambda value for input of translasso:
  lambda.const = 0
  ##initial
  if(ridge==F){
    if(family == "binomial"){
      # table(y1) %>% print
    }
    bini = glm(y1~x1,family=family)$coef}
  if(ridge==T){
    fit = cv.glmnet(x1, y1, alpha = 0, standardize = FALSE, lambda = lambda.grid, family = family)
    # Get the best lambda value
    best_lambda = fit$lambda.min
    ridge_fit = glmnet(x1, y1, alpha = 0, standardize = FALSE, lambda = best_lambda, family = family)
    # Get the coefficient vector
    bini = coef(ridge_fit) %>% as.vector()
    lambda.const = best_lambda
  }
  if(family=="binomial"){
    update = iteration.fun2(dat.list=dat.list,bini=bini,kk.list=1:K)
  }else{
    update = iteration.fun3(dat.list=dat.list,bini=bini,kk.list=1:K)
  }
  bnew.update=apply(update$b.k,1,mean) # updates need to be stopped here if using both binomial and Gaussian'
  betahat = bnew.update
  if(niter>1){
    for(ii in 1:(niter-1)){
      bnew=bnew.update
      if(family=="binomial"){
        #print("calling iteration.fun2")
        update.DCOS = iteration.fun2(dat.list=dat.list,bini=bnew,kk.list=1:K)}
      else{
        update.DCOS = iteration.fun3(dat.list=dat.list,bini=bnew,kk.list=1:K)}
      bnew.update = apply(update.DCOS$b.k,1,mean) 
    }
    Ahalf.DCOS= svd(-update.DCOS$Ahat); Ahalf.DCOS = Ahalf.DCOS$u%*%diag(sqrt(Ahalf.DCOS$d))%*%t(Ahalf.DCOS$v)
    res_list = Est.ALASSO.Approx.GLMNET(ynew=Ahalf.DCOS%*%bnew.update,xnew=Ahalf.DCOS,bini=bnew.update,N.adj=N,lambda.grid,family="gaussian")
    betahat = res_list$bhat.BIC
    lambda.BIC = res_list$lambda.BIC
    lambda.modBIC = res_list$ambda.modBIC
    lambda.const = lambda.BIC
  }
  res = list(coef = betahat, lam = lambda.const)
  return(res)
}

Est.ALASSO.Approx.GLMNET=function (ynew, xnew, bini, N.adj, lambda.grid, modBIC = T, N.inflate = N.adj, family) 
{
  w.b = 1/abs(bini)
  tmpfit = glmnet(x = xnew, y = ynew, family = family, 
                  penalty.factor = w.b, alpha = 1, lambda = lambda.grid, intercept = F)
  LL = apply((c(ynew) - predict(tmpfit, xnew, type = "response"))^2, 
             2, sum) * N.inflate
  if (modBIC) {
    BIC.lam = LL + min(N.adj^0.1, log(N.adj)) * tmpfit$df
    m.opt = which.min(BIC.lam)
    bhat.modBIC = tmpfit$beta[, m.opt]
    lamhat.modBIC = tmpfit$lambda[m.opt]
  }
  else {
    bhat.modBIC = lamhat.modBIC = NA
  }
  BIC.lam = LL + log(N.adj) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.BIC = tmpfit$beta[, m.opt]
  lamhat.BIC = tmpfit$lambda[m.opt]
  return(list(bhat.BIC = bhat.BIC, bhat.modBIC = bhat.modBIC, 
              lambda.BIC = lamhat.BIC, lambda.modBIC = lamhat.modBIC))
}

A.fun=function (bet, dat) 
{
  yy = dat[, 1]
  xx.vec = cbind(1, dat[, -1])
  -t(c(dg.logit(xx.vec %*% bet)) * xx.vec) %*% xx.vec/length(yy)
}
# for continuous var:
A.fun.2=function (bet, dat) 
{
  yy = dat[, 1]
  xx.vec = cbind(1, dat[, -1])
  2 * t(xx.vec) %*% xx.vec / length(yy)
}
U.fun=function (bet, dat) 
{
  yy = dat[, 1]
  xx.vec = cbind(1, dat[, -1])
  c(t(c(yy - g.logit(xx.vec %*% bet))) %*% xx.vec)/length(yy)
}
# for continuous var:
U.fun.2=function (bet, dat) 
{
  yy = dat[, 1]
  xx.vec = cbind(1, dat[, -1])
  -2 * t(xx.vec) %*% (yy - xx.vec %*% bet) / length(yy)
}
dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}

ROC.Est.FUN <- function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
{
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
  for(k in 1:pp)
  {
    yy = yy0; 
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }else{
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
  out
}

S.FUN <- function(yy,Yi,Di,yes.smooth=F)
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  }else{
    return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
}

sum.I <- function(yy,FUN,Yi,Vi=NULL){
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of descending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth)
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

# A self-defined matrix multiplication function that treats NA values as 0
matmul <- function(x, y) {ifelse(is.na(x), 0, x) %*% ifelse(is.na(y), 0, y)}
