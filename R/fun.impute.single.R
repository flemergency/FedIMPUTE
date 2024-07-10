# option 1: option = "rbinom"
imp.single.local.1 <- function(dat_original, var_name, family, myseed, ns, ns_option){
  maxit = 1
  # drop response y:
  dat_original = dat_original[, !(colnames(dat_original) == "y")]
  p = length(var_name)
  missing_rate = colMeans(is.na(dat_original[, var_name]))
  missing_rank = names(missing_rate)[order(missing_rate)] # rank var_name based on missing rate
  ### 1. use option = "rbinom"
  dat_temp = dat_original
  # first impute dat_temp by mean prediction:
  dat_temp = imp.mean(dat_temp, var_name, family)
  # store coefficients computed by the final round:
  coef_all = vector("list", length=p)
  for(i in 1:maxit){
    res = imp.mice.it1(dat_original = dat_original, dat_l = dat_temp, var_name = var_name, family = family, missing_rank = missing_rank, option = "rbinom", myseed=myseed, prev = NULL, ns, ns_option)
    dat_temp = res[[p]]$dat
    if(i==maxit){
      for(j in 1:p){
        coef_all[[j]] = res[[j]]$coef
      }
    }
  }
  return(list(dat = dat_temp, coef = coef_all))
}

# option2: prevalence
imp.single.local.2 <- function(dat_original, var_name, family, myseed, prevalence, ns, ns_option){
  maxit = 1
  # drop response y:
  dat_original = dat_original[, !(colnames(dat_original) == "y")]
  p = length(var_name)
  missing_rate = colMeans(is.na(dat_original[, var_name]))
  missing_rank = names(missing_rate)[order(missing_rate)] # rank var_name based on missing rate
  ### 2. use option = "prevalence"
  dat_temp = dat_original
  # first impute dat_temp by mean prediction:
  dat_temp = imp.mean(dat_temp, var_name, family)
  # store coefficients computed by the final round:
  coef_all = vector("list", length=p)
  for(i in 1:maxit){
    res = imp.mice.it1(dat_original = dat_original, dat_l = dat_temp, var_name = var_name, family = family, missing_rank = missing_rank, option = "prevalence", myseed=myseed, prev = prevalence, ns, ns_option)
    dat_temp = res[[p]]$dat
    if(i==maxit){
      for(j in 1:p){
        coef_all[[j]] = res[[j]]$coef
      }
    }
  }
  return(list(dat = dat_temp, coef = coef_all))
}

imp.single.local.all <- function(mylist, ns, ns_option){
  myseed = mylist$seed
  family = mylist$family
  dat_miss =mylist$dat_miss
  K = length(mylist$site_miss)
  col_miss = mylist$col_miss
######## option 1: rbinom
  # MNAR:
  imp_MNAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MNAR
  for(i in 1:K){
    res = imp.single.local.1(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], myseed, ns, ns_option)
    imp_MNAR[[i]]$dat = res$dat
    imp_MNAR[[i]]$coef = res$coef
  }
  # MAR:
  imp_MAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MAR
  for(i in 1:K){
    res = imp.single.local.1(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], myseed, ns, ns_option)
    imp_MAR[[i]]$dat = res$dat
    imp_MAR[[i]]$coef = res$coef
  }
  # MCAR:
  imp_MCAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MCAR
  for(i in 1:K){
    res = imp.single.local.1(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], myseed, ns, ns_option)
    imp_MCAR[[i]]$dat = res$dat
    imp_MCAR[[i]]$coef = res$coef
  }
  mylist$dat.imp.single.local.option1 = list(MNAR = imp_MNAR, MAR = imp_MAR, MCAR = imp_MCAR)

######## option 2: prevalence
  # first get prevalence: a vector of length(col_miss)
  prevalence = get_prevalence(mylist)
  # MNAR:
  imp_MNAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MNAR
  for(i in 1:K){
    res = imp.single.local.2(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], myseed, prevalence, ns, ns_option)
    imp_MNAR[[i]]$dat = res$dat
    imp_MNAR[[i]]$coef = res$coef
  }
  # MAR:
  imp_MAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MAR
  for(i in 1:K){
    res = imp.single.local.2(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], myseed, prevalence, ns, ns_option)
    imp_MAR[[i]]$dat = res$dat
    imp_MAR[[i]]$coef = res$coef
  }
  # MCAR:
  imp_MCAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MCAR
  for(i in 1:K){
    res = imp.single.local.2(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], myseed, prevalence, ns = F, ns_option=NULL)
    imp_MCAR[[i]]$dat = res$dat
    imp_MCAR[[i]]$coef = res$coef
  }
  mylist$dat.imp.single.local.option2 = list(MNAR = imp_MNAR, MAR = imp_MAR, MCAR = imp_MCAR)
  
  return(mylist)
}

########## pooled:
########## 
# pool coef:
get.imp.pool.single.coef <- function(mylist){
  K = mylist$K
  K0 = K - length(mylist$site_miss)
  var_miss = mylist$col_miss
  p = mylist$p
  sample_size = mylist$N
###### option 1: use rbinom()
  coef_meta = vector(mode = "list", length = 3)
  names(coef_meta) = c("MNAR", "MAR", "MCAR")
  coef_mean = vector(mode = "list", length = 3)
  names(coef_mean) = c("MNAR", "MAR", "MCAR")
  # loop through all var_miss:
  for(j in 1:length(var_miss)){
    # MNAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.single.local.option1$MNAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MNAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MNAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.single.local.option1$MAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MCAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.single.local.option1$MCAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MCAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MCAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
  }
  mylist$coef_single_pooled_option1 = list(coef_meta = coef_meta, coef_mean = coef_mean)

###### option 2: use prevalence 
  coef_meta = vector(mode = "list", length = 3)
  names(coef_meta) = c("MNAR", "MAR", "MCAR")
  coef_mean = vector(mode = "list", length = 3)
  names(coef_mean) = c("MNAR", "MAR", "MCAR")
  # loop through all var_miss:
  for(j in 1:length(var_miss)){
    # MNAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.single.local.option2$MNAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MNAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MNAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.single.local.option2$MAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MCAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.single.local.option2$MCAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MCAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MCAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
  }
  mylist$coef_single_pooled_option2 = list(coef_meta = coef_meta, coef_mean = coef_mean)  
  return(mylist)
}

# impute using pooled mice for all variables
# is a vector of length p
imp.single.pool.all<- function(mylist){
  K = mylist$K
  K0 = K - length(mylist$site_miss)
  var_name = mylist$col_miss
  family = mylist$family
####### option1: rbinom()
  dat_MNAR.meta = vector("list", length = K-K0)
  dat_MAR.meta = vector("list", length = K-K0)
  dat_MCAR.meta = vector("list", length = K-K0)
  dat_MNAR.mean = vector("list", length = K-K0)
  dat_MAR.mean = vector("list", length = K-K0)
  dat_MCAR.mean = vector("list", length = K-K0)
  # loop through each site with missingness:
  for(i in 1:(K-K0)){
    # MNAR:
    dat_original1 = mylist$dat_miss$MNAR[[i]] # original missingness
    dat_new1.meta = mylist$dat.imp.single.local.option1$MNAR[[i]]$dat
    dat_new1.mean = mylist$dat.imp.single.local.option1$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2.meta = mylist$dat.imp.single.local.option1$MAR[[i]]$dat
    dat_new2.mean = mylist$dat.imp.single.local.option1$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3.meta = mylist$dat.imp.single.local.option1$MCAR[[i]]$dat
    dat_new3.mean = mylist$dat.imp.single.local.option1$MCAR[[i]]$dat
      # loop through all missing variables:
      for(j in 1:length(var_name)){
      # MNAR:
      # meta
      coef = mylist$coef_single_pooled_option1$coef_meta$MNAR[[j]]
      dat_new1.meta = imp.mice.pool.1(dat_original1, dat_new1.meta, coef, var_name[j], family[j]) 
      # for the above imp.function, it does the same operation for both mice and single
      # mean:
      coef = mylist$coef_single_pooled_option1$coef_mean$MNAR[[j]]
      dat_new1.mean = imp.mice.pool.1(dat_original1, dat_new1.mean, coef, var_name[j], family[j])
      
      # MAR:
      # meta
      coef = mylist$coef_single_pooled_option1$coef_meta$MAR[[j]]
      dat_new2.meta = imp.mice.pool.1(dat_original2, dat_new2.meta, coef, var_name[j], family[j]) 
      # mean:
      coef = mylist$coef_single_pooled_option1$coef_mean$MAR[[j]]
      dat_new2.mean = imp.mice.pool.1(dat_original2, dat_new2.mean, coef, var_name[j], family[j])
      
      # MCAR:
      # meta
      coef = mylist$coef_single_pooled_option1$coef_meta$MCAR[[j]]
      dat_new3.meta = imp.mice.pool.1(dat_original3, dat_new3.meta, coef, var_name[j], family[j]) 
      # mean:
      coef = mylist$coef_single_pooled_option1$coef_mean$MCAR[[j]]
      dat_new3.mean = imp.mice.pool.1(dat_original3, dat_new3.mean, coef, var_name[j], family[j])
     
      }
    dat_MNAR.meta[[i]] = dat_new1.meta
    dat_MAR.meta[[i]] = dat_new2.meta
    dat_MCAR.meta[[i]] = dat_new3.meta
    
    dat_MNAR.mean[[i]] = dat_new1.mean
    dat_MAR.mean[[i]] = dat_new2.mean
    dat_MCAR.mean[[i]] = dat_new3.mean
  }
  mylist$dat.imp.single.pool.meta.option1 = list(MNAR = dat_MNAR.meta, MAR = dat_MAR.meta, MCAR = dat_MCAR.meta)
  mylist$dat.imp.single.pool.mean.option1 = list(MNAR = dat_MNAR.mean, MAR = dat_MAR.mean, MCAR = dat_MCAR.mean)

####### option2: prevalence
  dat_MNAR.meta = vector("list", length = K-K0)
  dat_MAR.meta = vector("list", length = K-K0)
  dat_MCAR.meta = vector("list", length = K-K0)
  dat_MNAR.mean = vector("list", length = K-K0)
  dat_MAR.mean = vector("list", length = K-K0)
  dat_MCAR.mean = vector("list", length = K-K0)
  # get prevalence
  prev = get_prevalence(mylist)
  # loop through each site with missingness:
  for(i in 1:(K-K0)){
    # MNAR:
    dat_original1 = mylist$dat_miss$MNAR[[i]] # original missingness
    dat_new1.meta = mylist$dat.imp.single.local.option2$MNAR[[i]]$dat
    dat_new1.mean = mylist$dat.imp.single.local.option2$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2.meta = mylist$dat.imp.single.local.option2$MAR[[i]]$dat
    dat_new2.mean = mylist$dat.imp.single.local.option2$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3.meta = mylist$dat.imp.single.local.option2$MCAR[[i]]$dat
    dat_new3.mean = mylist$dat.imp.single.local.option2$MCAR[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # MNAR:
      # meta
      coef = mylist$coef_single_pooled_option2$coef_meta$MNAR[[j]]
      dat_new1.meta = imp.mice.pool.2(dat_original1, dat_new1.meta, coef, var_name[j], family[j], prev[j]) 
      # mean:
      coef = mylist$coef_single_pooled_option2$coef_mean$MNAR[[j]]
      dat_new1.mean = imp.mice.pool.2(dat_original1, dat_new1.mean, coef, var_name[j], family[j], prev[j])
      
      # MAR:
      # meta
      coef = mylist$coef_single_pooled_option2$coef_meta$MAR[[j]]
      dat_new2.meta = imp.mice.pool.2(dat_original2, dat_new2.meta, coef, var_name[j], family[j], prev[j]) 
      # mean:
      coef = mylist$coef_single_pooled_option2$coef_mean$MAR[[j]]
      dat_new2.mean = imp.mice.pool.2(dat_original2, dat_new2.mean, coef, var_name[j], family[j], prev[j])
      
      # MCAR:
      # meta
      coef = mylist$coef_single_pooled_option2$coef_meta$MCAR[[j]]
      dat_new3.meta = imp.mice.pool.2(dat_original3, dat_new3.meta, coef, var_name[j], family[j], prev[j]) 
      # mean:
      coef = mylist$coef_single_pooled_option2$coef_mean$MCAR[[j]]
      dat_new3.mean = imp.mice.pool.2(dat_original3, dat_new3.mean, coef, var_name[j], family[j], prev[j])
      
    }
    dat_MNAR.meta[[i]] = dat_new1.meta
    dat_MAR.meta[[i]] = dat_new2.meta
    dat_MCAR.meta[[i]] = dat_new3.meta
    
    dat_MNAR.mean[[i]] = dat_new1.mean
    dat_MAR.mean[[i]] = dat_new2.mean
    dat_MCAR.mean[[i]] = dat_new3.mean
  }
  mylist$dat.imp.single.pool.meta.option2 = list(MNAR = dat_MNAR.meta, MAR = dat_MAR.meta, MCAR = dat_MCAR.meta)
  mylist$dat.imp.single.pool.mean.option2 = list(MNAR = dat_MNAR.mean, MAR = dat_MAR.mean, MCAR = dat_MCAR.mean)
  return(mylist)
}


# get coef estimated by DAC.
get.single.DAC <- function(mylist){
  # print("Doing DAC based on single imputation")
  var_name = mylist$col_miss
  K = mylist$K
  K0 = K- length(mylist$site_miss) # of sites with no missingness
  family = mylist$family
  niter = mylist$it_DAC
  ridge = mylist$ridge
  ###### option1: rbinom()
  # MNAR:
  coef.list1 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.single.local.option1$MNAR[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[i]] = cbind(y,X)
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = mylist$dat_original[[(K-K0)+i]]$dat
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[(K-K0)+i]] = cbind(y,X)
    }
    # get coef estimated by DAC:
    coef.list1[[j]] = DCOS.FUN(rev(dat.list), niter, ridge, NULL, family[j])
  }
  # print("finished MNAR")
  #MAR:
  coef.list2 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.single.local.option1$MAR[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[i]] = cbind(y,X)
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = mylist$dat_original[[(K-K0)+i]]$dat
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[(K-K0)+i]] = cbind(y,X)
    }
    # get coef estimated by DAC:
    coef.list2[[j]] = DCOS.FUN(rev(dat.list), niter, ridge, NULL, family[j])
  }
  #print("finished MAR")
  # MCAR
  coef.list3 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.single.local.option1$MCAR[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[i]] = cbind(y,X)
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = mylist$dat_original[[(K-K0)+i]]$dat
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[(K-K0)+i]] = cbind(y,X)
    }
    # get coef estimated by DAC:
    coef.list3[[j]] = DCOS.FUN(rev(dat.list), niter, ridge, lambda=NULL, family = family[j])
  }
  mylist$coef_DAC_single_option1 = list(MNAR = coef.list1, MAR = coef.list2, MCAR = coef.list3)
  
  
  ###### option2: prevalence
  # get prevalence"
  prev = get_prevalence(mylist)
  # MNAR:
  coef.list1 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.single.local.option2$MNAR[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[i]] = cbind(y,X)
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = mylist$dat_original[[(K-K0)+i]]$dat
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[(K-K0)+i]] = cbind(y,X)
    }
    # get coef estimated by DAC:
    # reverse list so that the no-missingness site is used for calculation of initial values
    coef.list1[[j]] = DCOS.FUN(rev(dat.list), niter, ridge, NULL, family[j])
  }
  #MAR:
  coef.list2 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.single.local.option2$MAR[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[i]] = cbind(y,X)
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = mylist$dat_original[[(K-K0)+i]]$dat
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[(K-K0)+i]] = cbind(y,X)
    }
    # get coef estimated by DAC:
    coef.list2[[j]] = DCOS.FUN(rev(dat.list), niter, ridge, NULL, family[j])
  }
  # MCAR
  coef.list3 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.single.local.option2$MCAR[[i]]$dat
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[i]] = cbind(y,X)
    }
    # now also join the dat of no-missingness sites:
    for(i in 1:K0){
      dat = mylist$dat_original[[(K-K0)+i]]$dat
      dat = dat[,-1] # dropping y response in needed as we're using dat_original
      dat2 = dat[,-which(colnames(dat) == var_name[j])] %>% as.matrix() # matrix for predictors of imputation model
      y_all = dat[,which(colnames(dat) == var_name[j])] %>% as.matrix() # set var_name as new response
      y = y_all[which(!is.na(y_all))]
      X = dat2[which(!is.na(y_all)),]
      dat.list[[(K-K0)+i]] = cbind(y,X)
    }
    # get coef estimated by DAC:
    coef.list3[[j]] = DCOS.FUN(rev(dat.list), niter, ridge, lambda=NULL, family = family[j])
  }
  mylist$coef_DAC_single_option2 = list(MNAR = coef.list1, MAR = coef.list2, MCAR = coef.list3)
  return(mylist)
}


imp.single.DAC.all <- function(mylist){
  K = mylist$K
  K0 = K - length(mylist$site_miss) # number of no-missingness sites
  var_name = mylist$col_miss
  family = mylist$family
  niter = mylist$it_DAC
  ####### option1: rbinom()
  dat_MNAR = vector("list", length = K-K0)
  dat_MAR = vector("list", length = K-K0)
  dat_MCAR = vector("list", length = K-K0)
  # loop through each site with missingness:
  for(i in 1:(K-K0)){
    # MNAR:
    dat_original1 = mylist$dat_miss$MNAR[[i]] # original missingness
    dat_new1 = mylist$dat.imp.single.local.option1$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2 = mylist$dat.imp.single.local.option1$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3 = mylist$dat.imp.single.local.option1$MCAR[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # MNAR:
      coef = mylist$coef_DAC_single_option1$MNAR[[j]]
      dat_new1 = imp.mice.pool.1(dat_original1, dat_new1, coef, var_name[j], family[j])
      # MAR:
      coef = mylist$coef_DAC_single_option1$MAR[[j]]
      dat_new2 = imp.mice.pool.1(dat_original2, dat_new2, coef, var_name[j], family[j])
      # MCAR:
      coef = mylist$coef_DAC_single_option1$MCAR[[j]]
      dat_new3 = imp.mice.pool.1(dat_original3, dat_new3, coef, var_name[j], family[j]) 
      
    }
    dat_MNAR[[i]] = dat_new1
    dat_MAR[[i]] = dat_new2
    dat_MCAR[[i]] = dat_new3
    
  }
  mylist$dat.imp.DAC.single.option1 = list(MNAR = dat_MNAR, MAR = dat_MAR, MCAR = dat_MCAR)
  ####### option2: prevalence
  dat_MNAR = vector("list", length = K-K0)
  dat_MAR = vector("list", length = K-K0)
  dat_MCAR = vector("list", length = K-K0)
  # get prevalence:
  prev = get_prevalence(mylist)
  # loop through each site:
  for(i in 1:(K-K0)){
    # MNAR:
    dat_original1 = mylist$dat_miss$MNAR[[i]] # original missingness
    dat_new1 = mylist$dat.imp.single.local.option2$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2 = mylist$dat.imp.single.local.option2$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3 = mylist$dat.imp.single.local.option2$MCAR[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # MNAR:
      coef = mylist$coef_DAC_single_option2$MNAR[[j]]
      dat_new1 = imp.mice.pool.2(dat_original1, dat_new1, coef, var_name[j], family[j], prev[j])
      # MAR:
      coef = mylist$coef_DAC_single_option2$MAR[[j]]
      dat_new2 = imp.mice.pool.2(dat_original2, dat_new2, coef, var_name[j], family[j], prev[j])
      # MCAR:
      coef = mylist$coef_DAC_single_option2$MCAR[[j]]
      dat_new3= imp.mice.pool.2(dat_original3, dat_new3, coef, var_name[j], family[j], prev[j])

    }
    dat_MNAR[[i]] = dat_new1
    dat_MAR[[i]] = dat_new2
    dat_MCAR[[i]] = dat_new3

  }
  mylist$dat.imp.DAC.single.option2 = list(MNAR = dat_MNAR, MAR = dat_MAR, MCAR = dat_MCAR)
  return(mylist)
}


