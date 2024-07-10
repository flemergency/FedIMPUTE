options(warn=-1)
# impute locally, using mean/median strategy 
imp.local.mean.all.realdata <- function(mylist, col_miss, family, site_miss){
  list = vector("list", length = length(site_miss))
  for(i in site_miss){
    dat_temp = mylist[[i]]
    list[[i]]= imp.mean.2(dat_temp, col_miss, family)
  }
  return(list)
}

# impute naive-federated, using mean/median from S3 (no-missingness site) for imputation
imp.naive.fed.all.realdata <- function(mylist, col_miss, family, site_miss){
  K = length(mylist)
  N = lapply(mylist, nrow) %>% unlist
  NN = N/sum(N)  # normalized weights
  K0 = length(site_miss)
  res0 = c()
  dat_temp = mylist[[K0+1]]
  for(i in 1:length(col_miss)){
    var_temp = dat_temp[[col_miss[i]]]
    res0_temp = ifelse(family[i]=="binomial",
                       as.numeric(names(sort(table(var_temp), decreasing = TRUE)[1])),
                       mean(var_temp))
    res0 = c(res0, res0_temp)
  }
  list = vector("list", length = length(site_miss))
  for(i in site_miss){
    dat_temp = mylist[[i]]
    dat_temp_new = imp.mean.3(dat_temp, col_miss, res0)
    list[[i]] = dat_temp_new
  }
  return(list)
}


imp.single.local.all.realdata <- function(dat_miss, col_miss, family, site_miss, myseed, ns, ns_option){
  K = length(site_miss)
  ######## option 1: rbinom
  imp1 = vector("list", length = K)
  dat_miss_temp = dat_miss
  for(i in 1:K){
    res = imp.single.local.1(dat_miss_temp[[i]], col_miss, family[1:length(col_miss)], myseed, ns, ns_option)
    imp1[[i]]$dat = res$dat
    imp1[[i]]$coef = res$coef
  }
  ######## option 2: prevalence
  # first get prevalence: a vector of length(col_miss)
  prevalence = get_prevalence_realdata(dat_miss, col_miss, family, site_miss)
  imp2 = vector("list", length = K)
  dat_miss_temp = dat_miss
  for(i in 1:K){
    res = imp.single.local.2(dat_miss_temp[[i]], col_miss, family[1:length(col_miss)], myseed, prevalence, ns, ns_option)
    imp2[[i]]$dat = res$dat
    imp2[[i]]$coef = res$coef
  }
  res = list(imp1, imp2)
  names(res) = c("option1", "option2")
  return(res)
}

imp.mice.local.all.realdata <- function(dat_miss, col_miss, family, site_miss, myseed, maxit = 5, ns, ns_option){
  K = length(site_miss)
  ######## option 1: rbinom
  imp1 = vector("list", length = K)
  dat_miss_temp = dat_miss
  for(i in 1:K){
    res = imp.mice.local.1(dat_miss_temp[[i]], col_miss, family[1:length(col_miss)], maxit, myseed, ns, ns_option)
    imp1[[i]]$dat = res$dat
    imp1[[i]]$coef = res$coef
  }
  ######## option 2: prevalence
  # first get prevalence: a vector of length(col_miss)
  prevalence = get_prevalence_realdata(dat_miss, col_miss, family, site_miss)
  imp2 = vector("list", length = K)
  dat_miss_temp = dat_miss
  for(i in 1:K){
    res = imp.mice.local.2(dat_miss_temp[[i]], col_miss, family[1:length(col_miss)], maxit, myseed, prevalence, ns, ns_option)
    imp2[[i]]$dat = res$dat
    imp2[[i]]$coef = res$coef
    imp2[[i]]$cutoff = res$cutoff  # the best cutoff for binary variable. 
  }
  res = list(imp1, imp2)
  names(res) = c("option1", "option2")
  return(res)
}

# get coef of non-missingness site for all variables
get.coef.nmiss.all.new.realdata <- function(dat_miss, col_miss, family, site_miss, ns, ns_option){
  K = length(dat_miss)
  var_miss = col_miss
  J = length(var_miss)
  dat = dat_miss
  coef_nmiss = vector("list", length = K)
  for(i in 1:K){
    coef = vector("list", length = J)
    dat_temp = dat[[i]]
    for(j in 1:J){
      coef_temp = get.coef.nmiss.single.realdata(var_miss[j], dat_temp, family[j], ns, ns_option)
      coef[[j]] = coef_temp
    }
    names(coef) = var_miss
    coef_nmiss[[i]] = coef
  }
  names(coef_nmiss) = paste0("S", c(1:K))
  return(coef_nmiss)
}

get.coef.nmiss.single.realdata <- function(var_name, dat, family, ns, ns_option){
  dat2 = dat[, -1] #drop response y
  dat2 = dat2 %>% na.omit()
  dat3 = subset(dat2, select = -which(names(dat2) == var_name)) %>% as.data.frame()
  y =  subset(dat2, select = which(names(dat2) == var_name)) %>% as.matrix()
  # add ns:
  # add ns for X:
  dat3_withns = dat3
  if(ns){
    for(i in 1:length(ns_option$ns_var)){
      # print(i)
      dat3_withns = cbind(dat3_withns, 
                             add_ns_onev(unlist(dat3[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
    }
  }
  # fit model with ns:
  res = lasso.fun(x = as.matrix(dat3_withns), y = as.vector(y), family = family)
  coef = res$coef
  return(coef)
}


########## pooled:
########## 
# pool coef:
get.imp.pool.coef.realdata <- function(dat_miss, dat_imp, col_miss, family, site_miss, myseed, ns, ns_option){
  K = length(dat_miss)
  site_all = c(1:K)
  var_miss = col_miss
  p = dat_miss[[1]] %>% ncol - 1
  sample_size = lapply(dat_miss, nrow) %>% unlist
  coef_non_missing = get.coef.nmiss.all.new.realdata(dat_miss, col_miss, family, site_miss, ns, ns_option)
  ###### option 1: use rbinom()
  # loop through all var_miss:
  for(j in 1:length(var_miss)){
    # print(j)
    v1 = c()
    for(i in site_miss){
      # print(i)
      # print(length(dat_imp$option1[[i]]$coef[[j]]))
      v1 = append(v1, dat_imp$option1[[i]]$coef[[j]])
    }
    # print(length(v1))
    for(i in site_all[-site_miss]){
      v1 = append(v1, coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    # print(length(v1))
    coef_meta1 = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean1 = get.coef.pool(v1, "mean", K, sample_size)
  }
  ###### option 2: use prevalence 
  # loop through all var_miss:
  for(j in 1:length(var_miss)){
    v1 = c()
    for(i in site_miss){
      v1 = append(v1, dat_imp$option2[[i]]$coef[[j]])
    }
    for(i in site_all[-site_miss]){
      v1 = append(v1,coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta2 = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean2 = get.coef.pool(v1, "mean", K, sample_size)
  }
  res = list(list(coef_meta1, coef_mean1), list(coef_meta2, coef_mean2))  
  names(res) = c("option1", "option2")
  return(res)
}

imp.pool.all.realdata<- function(dat_miss, dat_imp, col_miss, family, site_miss, myseed, ns, ns_option){
  K = length(dat_miss)
  K0 = K - length(site_miss)
  var_name = col_miss
  imp.coef = get.imp.pool.coef.realdata(dat_miss, dat_imp, col_miss, family, site_miss, myseed, ns, ns_option)
  # print(imp.coef)
  ####### option1: rbinom()
  dat.meta1 = vector("list", length = K-K0)
  dat.mean1 = vector("list", length = K-K0)
  # loop through each site with missingness:
  for(i in 1:(K-K0)){
    dat_original1 = dat_miss[[i]] 
    dat_new1.meta = dat_imp$option1[[i]]$dat
    dat_new1.mean = dat_imp$option1[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # meta
      coef = imp.coef$option1[[1]]
      dat_new1.meta = imp.mice.pool.1(dat_original1, dat_new1.meta, coef, var_name[j], family[j], ns, ns_option) 
      # for the above imp.function, it does the same operation for both mice and single
      # mean:
      coef = imp.coef$option1[[2]]
      dat_new1.mean = imp.mice.pool.1(dat_original1, dat_new1.mean, coef, var_name[j], family[j], ns, ns_option)
    }
    dat.meta1[[i]] = dat_new1.meta
    dat.mean1[[i]] = dat_new1.mean
  }
  ####### option2: prevalence
  # get prevalence
  prev = get_prevalence_realdata(dat_miss, col_miss, family, site_miss)
  dat.meta2 = vector("list", length = K-K0)
  dat.mean2 = vector("list", length = K-K0)
  # loop through each site with missingness:
  for(i in 1:(K-K0)){
    dat_original1 = dat_miss[[i]] 
    dat_new1.meta = dat_imp$option2[[i]]$dat
    dat_new1.mean = dat_imp$option2[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # meta
      coef = imp.coef$option1[[1]]
      dat_new1.meta = imp.mice.pool.1(dat_original1, dat_new1.meta, coef, var_name[j], family[j], ns, ns_option) 
      # for the above imp.function, it does the same operation for both mice and single
      # mean:
      coef = imp.coef$option1[[2]]
      dat_new1.mean = imp.mice.pool.1(dat_original1, dat_new1.mean, coef, var_name[j], family[j], ns, ns_option)
    }
    dat.meta2[[i]] = dat_new1.meta
    dat.mean2[[i]] = dat_new1.mean
  }
  
  res1 = list(dat.meta1, dat.mean1)
  names(res1) = c("meta", "mean")
  res2 = list(dat.meta2, dat.mean2)
  names(res2) = c("meta", "mean")
  res = list(res1, res2)
  names(res) = c("option1", "option2")
  
  return(res)
}

# for future development. need to give options to adjust the order of dat list
# get coef estimated by DAC.
# modified: return not only coef, but also lam.const for next TL
get.fl.realdata <- function(dat_miss, dat_imp, col_miss, family, site_miss, niter, ridge, ns, ns_option, DAC){
  var_name = col_miss
  K = length(dat_miss)
  K0 = K-length(site_miss)
  site_all = c(1:K)
  ###### option1: rbinom()
  coef.list1 = vector("list", length = length(var_name))
  lam.list1 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      # print(i)
      dat = dat_imp$option1[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),] %>% as.data.frame()
      # add ns for X:
      X_ns = X
      if(ns){
        for(m in 1:length(ns_option$ns_var)){
          # print(i)
          X_ns = cbind(X_ns, add_ns_onev(unlist(X[ns_option$ns_var[m]]), ns_option$ns_var[m], ns_option$df[m], ns_option$knots[[m]]))
        }
      }
      dat.list[[i]] = cbind(y,as.matrix(X_ns))
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = dat_miss[[(K-K0)+i]]
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),] %>% as.data.frame()
      # add ns for X:
      X_ns = X
      if(ns){
        for(m in 1:length(ns_option$ns_var)){
          # print(i)
          X_ns = cbind(X_ns, add_ns_onev(unlist(X[ns_option$ns_var[m]]), ns_option$ns_var[m], ns_option$df[m], ns_option$knots[[m]]))
        }
      }
      dat.list[[(K-K0)+i]] = cbind(y,as.matrix(X_ns))
    }
    if(DAC){
      # print("DAC results, option1:")
      # get coef estimated by DAC:
      res.temp = DCOS.FUN(rev(dat.list), niter, ridge, NULL, family[j])
    }
    else{
      # print("SHIR results, option1:")
      res.temp = SHIR_new_lm(dat.list, family[j])
    }
    # print(res.temp$coef)
    coef.list1[[j]] = res.temp$coef
    lam.list1[[j]] = res.temp$lam
  }
  ###### option2: prevalence
  # get prevalence"
  prevalence = get_prevalence_realdata(dat_miss, col_miss, family, site_miss)
  coef.list2 = vector("list", length = length(var_name))
  lam.list2 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      # print(i)
      dat = dat_imp$option2[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),] %>% as.data.frame()
      # add ns for X:
      X_ns = X
      if(ns){
        for(m in 1:length(ns_option$ns_var)){
          # print(i)
          X_ns = cbind(X_ns, add_ns_onev(unlist(X[ns_option$ns_var[m]]), ns_option$ns_var[m], ns_option$df[m], ns_option$knots[[m]]))
        }
      }
      dat.list[[i]] = cbind(y,as.matrix(X_ns))
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = dat_miss[[(K-K0)+i]]
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),] %>% as.data.frame()
      # add ns for X:
      X_ns = X
      if(ns){
        for(m in 1:length(ns_option$ns_var)){
          # print(i)
          X_ns = cbind(X_ns, add_ns_onev(unlist(X[ns_option$ns_var[m]]), ns_option$ns_var[m], ns_option$df[m], ns_option$knots[[m]]))
        }
      }
      dat.list[[(K-K0)+i]] = cbind(y,as.matrix(X_ns))
    }
    if(DAC){
      # print("DAC results, option2:")
      res.temp = DCOS.FUN(rev(dat.list), niter, ridge, NULL, family[j])
    }
    else{
      # print("SHIR results, option2:")
      res.temp = SHIR_new_lm(dat.list, family[j])
    }
    # print(res.temp$coef)
    coef.list2[[j]] = res.temp$coef
    lam.list2[[j]] = res.temp$lam
  }
  res = list(imp.coef = list(option1=coef.list1, option2=coef.list2), imp.lam = list(option1=lam.list1, option2=lam.list2))
  return(res)
}


imp.fl.all.realdata <- function(dat_miss, dat_imp, col_miss, family, site_miss, niter = 1, ridge = T, ns, ns_option, TL, DAC=T, penalty = F, lam.TL.target=F, n.s=NULL){
  var_name = col_miss
  K = length(dat_miss)
  K0 = K-length(site_miss)
  site_all = c(1:K)
  # fl:
  fed.res = get.fl.realdata(dat_miss, dat_imp, col_miss, family, site_miss, niter, ridge, ns, ns_option, DAC)
  
  imp.coef = fed.res$imp.coef
  imp.lam = fed.res$imp.lam
  
  ####### option1: rbinom()
  dat1_res = vector("list", length = K-K0)
  # loop through each site with missingness:
  for(i in 1:(K-K0)){
    dat_original = dat_miss[[i]] 
    dat_new1 = dat_imp$option1[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      coef = imp.coef$option1[[j]]
      # print("Original imputation model:")
      # print(coef)
      # add TL step to recalibrate coef:
      if(TL){
        # drop response y:
        dat_temp = dat_new1[, !(colnames(dat_new1) == "y")]
        # replace targeted variable back to missing status:
        dat_temp[[var_name[j]]]=dat_original[[var_name[j]]]
        dat2 = subset(dat_temp, select = -which(colnames(dat_temp) == var_name[j])) %>% as.matrix()
        y =  subset(dat_temp, select = which(colnames(dat_temp) == var_name[j])) %>% as.matrix()
        y.t= y[which(!is.na(y))]
        x.t = dat2[which(!is.na(y)),] %>% as.data.frame()
        x.t.new = x.t
        if(ns){
          for(i in 1:length(ns_option$ns_var)){
            x.t.new = cbind(x.t.new, add_ns_onev(unlist(x.t[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
          }
        }
        coef_TL = Trans.fun(coef[-1],as.matrix(x.t.new),y.t,n.s,ncol(x.t.new),imp.lam$option1[[j]],family[j],1-ridge, penalty, lam.TL.target)
        coef = coef_TL
      }
        dat_new1 = imp.mice.pool.1(dat_original, dat_new1, coef, var_name[j], family[j], ns, ns_option)
    }
    dat1_res[[i]] = dat_new1
  }
  ####### option2: prevalence
  dat2_res = vector("list", length = K-K0)
  # get prevalence:
  prev = get_prevalence_realdata(dat_miss, col_miss, family, site_miss)
  # loop through each site:
  for(i in 1:(K-K0)){
    dat_original = dat_miss[[i]] # original missingness
    dat_new1 = dat_imp$option2[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      coef = imp.coef$option2[[j]]
      # print("Original imputation model:")
      # print(coef)
      if(TL){
        # drop response y:
        dat_temp = dat_new1[, !(colnames(dat_new1) == "y")]
        # replace targeted variable back to missing status:
        dat_temp[[var_name[j]]]=dat_original[[var_name[j]]]
        dat2 = subset(dat_temp, select = -which(colnames(dat_temp) == var_name[j])) %>% as.matrix()
        y =  subset(dat_temp, select = which(colnames(dat_temp) == var_name[j])) %>% as.matrix()
        y.t= y[which(!is.na(y))]
        x.t = dat2[which(!is.na(y)),] %>% as.data.frame()
        x.t.new = x.t
        if(ns){
          for(i in 1:length(ns_option$ns_var)){
            x.t.new = cbind(x.t.new, add_ns_onev(unlist(x.t[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
          }
        }
        coef_TL = Trans.fun(coef[-1],as.matrix(x.t.new),y.t,n.s,ncol(x.t.new),imp.lam$option1[[j]],family[j],1-ridge, penalty, lam.TL.target)
        coef = coef_TL
      }
      dat_new1 = imp.mice.pool.2(dat_original, dat_new1, coef, var_name[j], family[j], prev[j], ns, ns_option)
    }
    dat2_res[[i]] = dat_new1
  }
  res = list(dat1_res, dat2_res)
  names(res) = c("option1", "option2")
  return(res)
  }

# Update: add ns for imputation model (except 1&2)
# 'penalty': for translasso
# lam.TL.target: option to fine tune lambda using target data
FedIMPUTE <- function(dat_miss, col_miss, family, site_miss, myseed, maxit_mice = 5, maxit_DAC = 1, ridge = T, ns = F, ns_option=NULL, penalty = F, lam.TL.target=F){
  # total sample size n.s (source sample size)
  n.s = sum(unlist(lapply(dat_miss, nrow)))
  # 1. local mean imputation
  dat_imp1 = imp.local.mean.all.realdata(dat_miss, col_miss, family, site_miss)
  # note: the $value is the mean of each local imputed variable, which may be used later
  print("1. Finished local mean imputation")
  # 2. naive federation
  dat_imp2 = imp.naive.fed.all.realdata(dat_miss, col_miss, family, site_miss)
  print("2. Finished naive imputation")
  # 3. local single imputation
  dat_imp3 = imp.single.local.all.realdata(dat_miss, col_miss, family, site_miss, myseed, ns, ns_option)
  print("3. Finished local single imputation")
  # note: 'option1' is the result by generating 0,1 using rbinom(); 'option2' is result by using prevalence; same for all following; no difference for continuous missing vars.
  # 4. local mice imputation
  dat_imp4 = imp.mice.local.all.realdata(dat_miss, col_miss, family, site_miss, myseed, maxit_mice, ns, ns_option)
  print("4. Finished local mice imputation")
  # 5. pooled single
  dat_imp5 = imp.pool.all.realdata(dat_miss, dat_imp3, col_miss, family, site_miss, myseed, ns, ns_option)
  print("5. Finished FedIMPUTE-Pooled single imputation")
  # 6. pooled mice
  dat_imp6 = imp.pool.all.realdata(dat_miss, dat_imp4, col_miss, family, site_miss, myseed, ns, ns_option)
  print("6. Finished FedIMPUTE-Pooled mice imputation")
  # 7. DAC-single
  dat_imp7 = imp.fl.all.realdata(dat_miss, dat_imp3, col_miss, family, site_miss, maxit_DAC, ridge, ns, ns_option, TL=F, DAC=T, penalty)
  print("7. Finished FedIMPUTE-DAC single imputation")
  # 8. DAC-mice
  dat_imp8 = imp.fl.all.realdata(dat_miss, dat_imp4, col_miss, family, site_miss, maxit_DAC, ridge, ns, ns_option, TL=F, DAC=T, penalty)
  print("8. Finished FedIMPUTE-DAC mice imputation")
  # 9. DAC-TL-single
  dat_imp9 = imp.fl.all.realdata(dat_miss, dat_imp3, col_miss, family, site_miss, maxit_DAC, ridge, ns, ns_option, TL=T, DAC=T, penalty, lam.TL.target, n.s)
  print("9. Finished FedIMPUTE-DAC single imputation with TL calibration")
  # 10. DAC-TL-MICE
  dat_imp10 = imp.fl.all.realdata(dat_miss, dat_imp4, col_miss, family, site_miss, maxit_DAC, ridge, ns, ns_option, TL=T, DAC=T, penalty, lam.TL.target, n.s)
  print("10. Finished FedIMPUTE-DAC mice imputation with TL calibration")
  
  res = list(dat_imp1, dat_imp2, dat_imp3, dat_imp4,
             dat_imp5, dat_imp6, dat_imp7, dat_imp8,
             dat_imp9, dat_imp10)
  names(res) = c("local mean", "naive fed", "local single", "local mice", 
                 "FedIMPUTE-pooled single", "FedIMPUTE-pooled mice", "FedIMPUTE-DAC single", "FedIMPUTE-DAC mice", 
                 "FedIMPUTE-DAC-TL single", "FedIMPUTE-DAC-TL mice")
  res
}

