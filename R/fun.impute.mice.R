library(splines)

imp.mean <- function(dat, var_name, family){
  for(i in 1:length(var_name)){
    var = dat[[var_name[i]]]
    if(family[i]=="binomial"){
      # print(sort(table(var)))
      freq_x <- as.numeric(names(sort(table(var), decreasing = TRUE)[1]))
      var[is.na(var)] = freq_x
    }
    else{
      var[is.na(var)] <- mean(var, na.rm = TRUE)
    }
    dat[[var_name[i]]] = var
  }
  return(dat)
}

# also return the freq/mean value
imp.mean.2 <- function(dat, var_name, family){
  vec = c()
  for(i in 1:length(var_name)){
    var = dat[[var_name[i]]]
    value = 0
    if(family[i]=="binomial"){
      freq_x <- as.numeric(names(sort(table(var), decreasing = TRUE)[1]))
      value = freq_x
      var[is.na(var)] = value
    }
    else{
      value = mean(var, na.rm = TRUE)
      var[is.na(var)] <- value
    }
    dat[[var_name[i]]] = var
    vec = append(vec, value)
  }
  names(vec) = var_name
  return(list(dat=dat, value = vec))
}

# use meta/mean pooled value to impute
imp.mean.3 <- function(dat, var_name, value){
  for(i in 1:length(var_name)){
    var = dat[[var_name[i]]]
    var[is.na(var)] = value[i]
    dat[[var_name[i]]] = var
  }
  return(dat)
}

# add ns for one variable:
add_ns_onev <- function(x, var_name, df, knots){
  # print(knots)
  res_ns = ns(x, df=df, knots=knots)
  res = as.data.frame(res_ns)
  colnames(res) = paste0(paste0(var_name, "ns"), c(1:ncol(res)))
  res
}

# impute a missing variable using the rest of the variables
# return a list with imputed data matrix and coef of imputation
# now instead of using cutoff to generate 0 or 1, we use
# rbinom()
imp.single.2 <- function(dat, var_name, family, myseed, ns, ns_option){
  set.seed(myseed)
  names = colnames(dat)
  dat2 = subset(dat, select = -which(names(dat) == var_name)) %>% as.matrix()
  y =  subset(dat, select = which(names(dat) == var_name)) %>% as.matrix()
  
  # split into training and testing data
  y_train = y[which(!is.na(y))]
  x_train = dat2[which(!is.na(y)),] %>% as.data.frame()
  x_test = dat2[which(is.na(y)),] 
  # X_test = cbind(rep(1, nrow(x_test)), x_test) %>% as.data.frame()
  # automatically match dimensions:
  if(!(ncol(x_test) == length(coef))){
    x_test = cbind(rep(1, nrow(x_test)), x_test) %>% as.data.frame()
  }
  
  # add ns for X:
  x_train_withns = x_train
  x_test_withns = x_test 
  if(ns){
    for(i in 1:length(ns_option$ns_var)){
      # print(i)
      x_train_withns = cbind(x_train_withns, 
                 add_ns_onev(unlist(x_train[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
      x_test_withns = cbind(x_test_withns, 
                       add_ns_onev(unlist(x_test[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
    }
  }
  
  # fit model:
  res = lasso.fun(x = x_train_withns, y = y_train, family = family)
  coef = res$coef
  pred = matmul(as.matrix(x_test_withns), coef)
  if (family == "binomial") {
    pred = tryCatch(
      {
        sapply(g.logit(pred), function(p) rbinom(1, 1, p))
      },
      warning = function(w) {
        message("Warning occurred during binomial prediction:", w)
        return(NA)
      }
    )
  }
  y[which(is.na(y))] = pred # replace NAs with predicted values
  y = y %>% as.data.frame()
  colnames(y) = var_name
  dat_imputed = cbind(dat2, y)
  dat_imputed = dat_imputed[, names]
  
  return(list(dat = dat_imputed, coef = coef))
}

# impute a missing variable using the rest of the variables
# return a list with imputed data matrix and coef of imputation
# now instead of using cutoff to generate 0 or 1, we use
# prevalence of this variable in the no-missingness site
imp.single.3 <- function(dat, var_name, family, myseed, prev, ns, ns_option){
  set.seed(myseed)
  names = colnames(dat)
  dat2 = subset(dat, select = -which(names(dat) == var_name)) %>% as.matrix()
  y =  subset(dat, select = which(names(dat) == var_name)) %>% as.matrix()
  
  # split into training and testing data
  y_train = y[which(!is.na(y))]
  x_train = dat2[which(!is.na(y)),] %>% as.data.frame()
  x_test = dat2[which(is.na(y)),]
  # automatically match dimensions:
  if(!(ncol(x_test) == length(coef))){
    x_test = cbind(rep(1, nrow(x_test)), x_test) %>% as.data.frame()
  }
  # add ns for X:
  x_train_withns = x_train
  x_test_withns = x_test 
  if(ns){
    for(i in 1:length(ns_option$ns_var)){
      x_train_withns = cbind(x_train_withns, 
                             add_ns_onev(unlist(x_train[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
      x_test_withns = cbind(x_test_withns, 
                            add_ns_onev(unlist(x_test[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
    }
  }
  
  res = lasso.fun(x = x_train_withns, y = y_train, family = family)
  coef = res$coef
  pred = matmul(as.matrix(x_test_withns), coef)
  if(family=="binomial"){
    cutoff = quantile(pred, probs = prev)
    pred = ifelse(pred>cutoff, 1, 0)
  }
  y[which(is.na(y))] = pred # replace NAs with predicted values
  
  names(coef) = c("Intercept", colnames(dat2))
  y = y %>% as.data.frame()
  colnames(y) = var_name
  dat_imputed = cbind(dat2, y)
  dat_imputed = dat_imputed[, names]
  return(list(dat = dat_imputed, coef = coef))
}

# dat_original: original data matrix
# dat_l: data matrix imputed at last round of iteration
# both var_name and family are vectors of length p; note that the two vectors should match
# missing_rank: a vector of length p which rank variables by missing ratio, from lowest to highest
# only one round of iterations
# return a list of length(var_name), say p
# then list[[p]] contains the final imputed data, and the coefficients of the last imputation model.
# To access the p imputation model for all p variables in this iteration,
# access list[[1]]$coef - list[[p]]$coef accordingly
imp.mice.it1<- function(dat_original, dat_l, var_name, family, missing_rank = NULL, option = "rbinom", myseed, prev=NULL, ns, ns_option){
  if(length(var_name)!=length(family)){
    print("error! length of var_name and family does not match")
    return()
  }
  dat_temp = dat_l
  res <- vector(mode = "list", length = length(var_name))
  # print("---------")
  # print(missing_rank)
  for(i in 1:length(missing_rank)){
    # print(missing_rank[i])
    # first change targeted variable back to missing status:
    dat_temp[[missing_rank[i]]] = dat_original[[missing_rank[i]]]
    # impute:
    if(option == "cutoff"){
      res_temp = imp.single(dat_temp, var_name = missing_rank[i],family = family[which(var_name == missing_rank[i])])
    }
    else if(option == "rbinom"){
      res_temp = imp.single.2(dat_temp, var_name = missing_rank[i],family = family[which(var_name == missing_rank[i])], myseed = myseed, ns, ns_option)
    }
    else{
      res_temp = imp.single.3(dat_temp, var_name = missing_rank[i],family = family[which(var_name == missing_rank[i])], myseed = myseed, prev = prev[which(var_name == missing_rank[i])], ns, ns_option)
    }
    dat_temp = res_temp$dat
    res[[i]]$dat = res_temp$dat
    res[[i]]$coef = res_temp$coef
    cut_off_temp = res_temp$cutoff
    res[[i]]$cutoff = cut_off_temp
  }
  return(res)
}

# option 1: option = "rbinom"
imp.mice.local.1 <- function(dat_original, var_name, family, maxit, myseed, ns, ns_option){
  # drop response y:
  dat_original = dat_original[, !(colnames(dat_original) == "y")]
  p = length(var_name)
  # print(colnames(dat_original))
  # print(var_name)
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
imp.mice.local.2 <- function(dat_original, var_name, family, maxit, myseed, prevalence, ns, ns_option){
  # drop response y:
  dat_original = dat_original[, !(colnames(dat_original) == "y")]
  p = length(var_name)
  # print(colnames(dat_original))
  # print(var_name)
  missing_rate = colMeans(is.na(dat_original[, var_name]))
  missing_rank = names(missing_rate)[order(missing_rate)] # rank var_name based on missing rate
  ### 2. use option = "prevalance"
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

imp.mice.local.all <- function(mylist, ns, ns_option){
  myseed = mylist$seed
  family = mylist$family
  maxit = mylist$it_mice
  dat_miss =mylist$dat_miss
  K = length(mylist$site_miss)
  col_miss = mylist$col_miss
######## option 1: rbinom
  # MNAR:
  imp_MNAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MNAR
  for(i in 1:K){
    res = imp.mice.local.1(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], maxit, myseed, ns, ns_option)
    imp_MNAR[[i]]$dat = res$dat
    imp_MNAR[[i]]$coef = res$coef
  }
  # MAR:
  imp_MAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MAR
  for(i in 1:K){
    res = imp.mice.local.1(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], maxit, myseed, ns, ns_option)
    imp_MAR[[i]]$dat = res$dat
    imp_MAR[[i]]$coef = res$coef
  }
  # MCAR:
  imp_MCAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MCAR
  for(i in 1:K){
    res = imp.mice.local.1(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], maxit, myseed, ns, ns_option)
    imp_MCAR[[i]]$dat = res$dat
    imp_MCAR[[i]]$coef = res$coef
  }
  mylist$dat.imp.mice.local.option1 = list(MNAR = imp_MNAR, MAR = imp_MAR, MCAR = imp_MCAR)

######## option 2: prevalence
  # first get prevalence: a vector of length(col_miss)
  prevalence = get_prevalence(mylist)
  # MNAR:
  imp_MNAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MNAR
  for(i in 1:K){
    res = imp.mice.local.2(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], maxit, myseed, prevalence, ns, ns_option)
    imp_MNAR[[i]]$dat = res$dat
    imp_MNAR[[i]]$coef = res$coef
    imp_MNAR[[i]]$cutoff = res$cutoff  # the best cutoff for binary variable.
  }
  # MAR:
  imp_MAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MAR
  for(i in 1:K){
    res = imp.mice.local.2(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], maxit, myseed, prevalence, ns, ns_option)
    imp_MAR[[i]]$dat = res$dat
    imp_MAR[[i]]$coef = res$coef
    imp_MAR[[i]]$cutoff = res$cutoff  # the best cutoff for binary variable.
  }
  # MCAR:
  imp_MCAR = vector("list", length = K)
  dat_miss_temp = dat_miss$MCAR
  for(i in 1:K){
    res = imp.mice.local.2(dat_miss_temp[[i]], col_miss, (mylist$family)[1:length(col_miss)], maxit, myseed, prevalence)
    imp_MCAR[[i]]$dat = res$dat
    imp_MCAR[[i]]$coef = res$coef
    imp_MCAR[[i]]$cutoff = res$cutoff  # the best cutoff for binary variable.
  }
  mylist$dat.imp.mice.local.option2 = list(MNAR = imp_MNAR, MAR = imp_MAR, MCAR = imp_MCAR)
  
  return(mylist)
}

########## pooled:
########## 
# pool coef:
get.imp.pool.mice.coef <- function(mylist){
  K = mylist$K
  K0 = K - length(mylist$site_miss) # K0: # of no-missingness site
  var_miss = mylist$col_miss
  p = mylist$p
  sample_size = mylist$N
  sample_size = sample_size[1:K]
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
      v1 = append(v1, mylist$dat.imp.mice.local.option1$MNAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MNAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MNAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.mice.local.option1$MAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MCAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.mice.local.option1$MCAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MCAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MCAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
  }
  mylist$coef_mice_pooled_option1 = list(coef_meta = coef_meta, coef_mean = coef_mean)

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
      v1 = append(v1, mylist$dat.imp.mice.local.option2$MNAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MNAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MNAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.mice.local.option2$MAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
    # MCAR:
    v1 = c()
    for(i in 1:(K-K0)){
      v1 = append(v1, mylist$dat.imp.mice.local.option2$MCAR[[i]]$coef[[j]])
    }
    for(i in 1:K0){
      v1 = append(v1, mylist$coef_non_missing[[i]][[j]])
    } # append v1 with coef of non-missingness site
    coef_meta$MCAR[[j]] = get.coef.pool(v1, "meta", K, sample_size)
    coef_mean$MCAR[[j]] = get.coef.pool(v1, "mean", K, sample_size)
  }
  mylist$coef_mice_pooled_option2 = list(coef_meta = coef_meta, coef_mean = coef_mean)  
  return(mylist)
}

get.coef.pool <- function(v1, method = "meta", K, sample_size){
  # print(length(v1))
  Coef = matrix(v1, byrow = F, ncol = K)
  if(method=="meta"){
    #print("Pooling by sample size weighted mean")
    coef = apply(Coef, 1, FUN = function(x){
      weighted.mean(x, w = sample_size/sum(sample_size))
    })
  }
  else{
    #print("Pooling by mean")
    coef = rowMeans(Coef)
  }
  return(coef)
}

# pool cutoff:
get.imp.pool.cutoff <- function(mylist){
  K = mylist$K
  var_miss = mylist$col_miss
  p = mylist$p
  sample_size = mylist$N
  cutoff_meta = vector(mode = "list", length = K)
  names(cutoff_meta) = c("MNAR", "MAR", "MCAR")
  cutoff_mean = vector(mode = "list", length = K)
  names(cutoff_mean) = c("MNAR", "MAR", "MCAR")
  # loop through all var_miss:
  for(j in 1:length(var_miss)){
    # MNAR:
    v1 = c()
    for(i in 1:K){
      v1 = append(v1, mylist$dat.imp.mice.local$MNAR[[i]]$cutoff[[j]])
    }
    cutoff_meta$MNAR[[j]] = get.cutoff.pool(v1, "meta", K, sample_size)
    cutoff_mean$MNAR[[j]] = get.cutoff.pool(v1, "mean", K, sample_size)
    # MAR:
    v1 = c()
    for(i in 1:K){
      v1 = append(v1, mylist$dat.imp.mice.local$MAR[[i]]$cutoff[[j]])
    }
    cutoff_meta$MAR[[j]] = get.cutoff.pool(v1, "meta", K, sample_size)
    cutoff_mean$MAR[[j]] = get.cutoff.pool(v1, "mean", K, sample_size)
    # MCAR:
    v1 = c()
    for(i in 1:K){
      v1 = append(v1, mylist$dat.imp.mice.local$MCAR[[i]]$cutoff[[j]])
    }
    cutoff_meta$MCAR[[j]] = get.cutoff.pool(v1, "meta", K, sample_size)
    cutoff_mean$MCAR[[j]] = get.cutoff.pool(v1, "mean", K, sample_size)
  }
  mylist$cutoff_mice_pooled = list(cutoff_meta = cutoff_meta, cutoff_mean = cutoff_mean)
  return(mylist)
}
get.cutoff.pool <- function(v1, method = "meta", K, sample_size){
  cutoff = matrix(v1, byrow = F, ncol = K)
  if(method=="meta"){
    res = apply(cutoff, 1, FUN = function(x){
      weighted.mean(x, w = sample_size/sum(sample_size))
    })
  }
  else{
    #print("Pooling by mean")
    res = rowMeans(cutoff)
  }
  return(res)
}

# option1: rbinom()
imp.mice.pool.1 <- function(dat_original, dat_l, coef, var_name, family, ns, ns_option){
  dat_temp = dat_l
  # drop response y:
  dat_temp = dat_temp[, !(colnames(dat_temp) == "y")]
  names = colnames(dat_temp)
  # replace targeted variable back to missing status:
  dat_temp[[var_name]] = dat_original[[var_name]]
  
  dat2 = subset(dat_temp, select = -which(colnames(dat_temp) == var_name)) %>% as.matrix()
  y =  subset(dat_temp, select = which(colnames(dat_temp) == var_name)) %>% as.matrix()
  
  y_train = y[which(!is.na(y))]
  x_train = dat2[which(!is.na(y)),]
  x_test = dat2[which(is.na(y)),]
  # automatically match dimensions:
  # this is not sufficient for added ns
  if(!(ncol(x_test) == length(coef))){
    x_test = cbind(rep(1, nrow(x_test)), x_test) %>% as.data.frame()
  }
  # add ns for X:
  x_test_withns = x_test 
  if(ns){
    for(i in 1:length(ns_option$ns_var)){
      x_test_withns = cbind(x_test_withns, 
                            add_ns_onev(unlist(x_test[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
    }
  }
  y_pred = as.matrix(x_test_withns) %*% coef
  
  if (family == "binomial") {
    pred = tryCatch(
      {
        sapply(g.logit(y_pred), function(p) rbinom(1, 1, p))
      },
      warning = function(w) {
        message("Warning occurred during binomial prediction:", w)
        return(NA)
      }
    )
  }
  y[which(is.na(y))] = y_pred # replace NAs with predicted values
  
  y = y %>% as.data.frame()
  colnames(y) = var_name
  dat_imputed = cbind(dat2, y)
  dat_imputed = dat_imputed[, names]
  
  return(dat_imputed)
}

# option2: prevalence
imp.mice.pool.2 <- function(dat_original, dat_l, coef, var_name, family, prev, ns, ns_option){
  dat_temp = dat_l
  # drop response y:
  dat_temp = dat_temp[, !(colnames(dat_temp) == "y")]
  names = colnames(dat_temp)
  # replace targeted variable back to missing status:
  dat_temp[[var_name]] = dat_original[[var_name]]
  
  dat2 = subset(dat_temp, select = -which(colnames(dat_temp) == var_name)) %>% as.matrix()
  y =  subset(dat_temp, select = which(colnames(dat_temp) == var_name)) %>% as.matrix()
  
  # split into training and testing data
  y_train = y[which(!is.na(y))]
  x_train = dat2[which(!is.na(y)),]
  x_test = dat2[which(is.na(y)),]
  # automatically match dimensions:
  if(!(ncol(x_test) == length(coef))){
    x_test = cbind(rep(1, nrow(x_test)), x_test) %>% as.data.frame()
  }
  # add ns for X:
  x_test_withns = x_test 
  if(ns){
    for(i in 1:length(ns_option$ns_var)){
      x_test_withns = cbind(x_test_withns, 
                            add_ns_onev(unlist(x_test[ns_option$ns_var[i]]), ns_option$ns_var[i], ns_option$df[i], ns_option$knots[[i]]))
    }
  }
  y_pred = as.matrix(x_test_withns) %*% coef
  if(family=="binomial"){
    cutoff = quantile(y_pred, probs = prev)
    y_pred = ifelse(y_pred>cutoff, 1, 0)
  }
  y[which(is.na(y))] = y_pred # replace NAs with predicted values
  y = y %>% as.data.frame()
  colnames(y) = var_name
  dat_imputed = cbind(dat2, y)
  dat_imputed = dat_imputed[, names]
  
  return(dat_imputed)
}

# impute using pooled mice for all variables
# is a vector of length p
imp.mice.pool.all<- function(mylist){
  K = mylist$K
  K0 = length(mylist$site_miss)
  var_name = mylist$col_miss
  family = mylist$family
####### option1: rbinom()
  dat_MNAR.meta = vector("list", length = K0)
  dat_MAR.meta = vector("list", length = K0)
  dat_MCAR.meta = vector("list", length = K0)
  dat_MNAR.mean = vector("list", length = K0)
  dat_MAR.mean = vector("list", length = K0)
  dat_MCAR.mean = vector("list", length = K0)
  # loop through each site with missingness:
  for(i in 1:K0){
    # MNAR:
    dat_original1 = mylist$dat_miss$MNAR[[i]] # original missingness
    dat_new1.meta = mylist$dat.imp.mice.local.option1$MNAR[[i]]$dat
    dat_new1.mean = mylist$dat.imp.mice.local.option1$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2.meta = mylist$dat.imp.mice.local.option1$MAR[[i]]$dat
    dat_new2.mean = mylist$dat.imp.mice.local.option1$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3.meta = mylist$dat.imp.mice.local.option1$MCAR[[i]]$dat
    dat_new3.mean = mylist$dat.imp.mice.local.option1$MCAR[[i]]$dat
      # loop through all missing variables:
      for(j in 1:length(var_name)){
      # MNAR:
      # meta
      coef = mylist$coef_mice_pooled_option1$coef_meta$MNAR[[j]]
      dat_new1.meta = imp.mice.pool.1(dat_original1, dat_new1.meta, coef, var_name[j], family[j]) 
      # mean:
      coef = mylist$coef_mice_pooled_option1$coef_mean$MNAR[[j]]
      dat_new1.mean = imp.mice.pool.1(dat_original1, dat_new1.mean, coef, var_name[j], family[j])
      
      # MAR:
      # meta
      coef = mylist$coef_mice_pooled_option1$coef_meta$MAR[[j]]
      dat_new2.meta = imp.mice.pool.1(dat_original2, dat_new2.meta, coef, var_name[j], family[j]) 
      # mean:
      coef = mylist$coef_mice_pooled_option1$coef_mean$MAR[[j]]
      dat_new2.mean = imp.mice.pool.1(dat_original2, dat_new2.mean, coef, var_name[j], family[j])
      
      # MCAR:
      # meta
      coef = mylist$coef_mice_pooled_option1$coef_meta$MCAR[[j]]
      dat_new3.meta = imp.mice.pool.1(dat_original3, dat_new3.meta, coef, var_name[j], family[j]) 
      # mean:
      coef = mylist$coef_mice_pooled_option1$coef_mean$MCAR[[j]]
      dat_new3.mean = imp.mice.pool.1(dat_original3, dat_new3.mean, coef, var_name[j], family[j])
     
      }
    dat_MNAR.meta[[i]] = dat_new1.meta
    dat_MAR.meta[[i]] = dat_new2.meta
    dat_MCAR.meta[[i]] = dat_new3.meta
    
    dat_MNAR.mean[[i]] = dat_new1.mean
    dat_MAR.mean[[i]] = dat_new2.mean
    dat_MCAR.mean[[i]] = dat_new3.mean
  }
  mylist$dat.imp.mice.pool.meta.option1 = list(MNAR = dat_MNAR.meta, MAR = dat_MAR.meta, MCAR = dat_MCAR.meta)
  mylist$dat.imp.mice.pool.mean.option1 = list(MNAR = dat_MNAR.mean, MAR = dat_MAR.mean, MCAR = dat_MCAR.mean)

####### option2: prevalence
  dat_MNAR.meta = vector("list", length = K0)
  dat_MAR.meta = vector("list", length = K0)
  dat_MCAR.meta = vector("list", length = K0)
  dat_MNAR.mean = vector("list", length = K0)
  dat_MAR.mean = vector("list", length = K0)
  dat_MCAR.mean = vector("list", length = K0)
  # get prevalence
  prev = get_prevalence(mylist)
  # loop through each site with missingness:
  for(i in 1:K0){
    # MNAR:
    dat_original1 = mylist$dat_miss$MNAR[[i]] # original missingness
    dat_new1.meta = mylist$dat.imp.mice.local.option2$MNAR[[i]]$dat
    dat_new1.mean = mylist$dat.imp.mice.local.option2$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2.meta = mylist$dat.imp.mice.local.option2$MAR[[i]]$dat
    dat_new2.mean = mylist$dat.imp.mice.local.option2$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3.meta = mylist$dat.imp.mice.local.option2$MCAR[[i]]$dat
    dat_new3.mean = mylist$dat.imp.mice.local.option2$MCAR[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # MNAR:
      # meta
      coef = mylist$coef_mice_pooled_option2$coef_meta$MNAR[[j]]
      dat_new1.meta = imp.mice.pool.2(dat_original1, dat_new1.meta, coef, var_name[j], family[j], prev[j]) 
      # mean:
      coef = mylist$coef_mice_pooled_option2$coef_mean$MNAR[[j]]
      dat_new1.mean = imp.mice.pool.2(dat_original1, dat_new1.mean, coef, var_name[j], family[j], prev[j])
      
      # MAR:
      # meta
      coef = mylist$coef_mice_pooled_option2$coef_meta$MAR[[j]]
      dat_new2.meta = imp.mice.pool.2(dat_original2, dat_new2.meta, coef, var_name[j], family[j], prev[j]) 
      # mean:
      coef = mylist$coef_mice_pooled_option2$coef_mean$MAR[[j]]
      dat_new2.mean = imp.mice.pool.2(dat_original2, dat_new2.mean, coef, var_name[j], family[j], prev[j])
      
      # MCAR:
      # meta
      coef = mylist$coef_mice_pooled_option2$coef_meta$MCAR[[j]]
      dat_new3.meta = imp.mice.pool.2(dat_original3, dat_new3.meta, coef, var_name[j], family[j], prev[j]) 
      # mean:
      coef = mylist$coef_mice_pooled_option2$coef_mean$MCAR[[j]]
      dat_new3.mean = imp.mice.pool.2(dat_original3, dat_new3.mean, coef, var_name[j], family[j], prev[j])
      
    }
    dat_MNAR.meta[[i]] = dat_new1.meta
    dat_MAR.meta[[i]] = dat_new2.meta
    dat_MCAR.meta[[i]] = dat_new3.meta
    
    dat_MNAR.mean[[i]] = dat_new1.mean
    dat_MAR.mean[[i]] = dat_new2.mean
    dat_MCAR.mean[[i]] = dat_new3.mean
  }
  mylist$dat.imp.mice.pool.meta.option2 = list(MNAR = dat_MNAR.meta, MAR = dat_MAR.meta, MCAR = dat_MCAR.meta)
  mylist$dat.imp.mice.pool.mean.option2 = list(MNAR = dat_MNAR.mean, MAR = dat_MAR.mean, MCAR = dat_MCAR.mean)
  return(mylist)
}


# get coef estimated by DAC.
get.mice.DAC <- function(mylist){
  # print("Doing DAC based on MICE imputation")
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
      dat = mylist$dat.imp.mice.local.option1$MNAR[[i]]$dat
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
      dat = mylist$dat.imp.mice.local.option1$MAR[[i]]$dat
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
    coef.list2[[j]] = DCOS.FUN(dat.list, niter, ridge, NULL, family[j])
  }
  # MCAR
  coef.list3 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.mice.local.option1$MCAR[[i]]$dat
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
  # print("finished MCAR")
  mylist$coef_DAC_mice_option1 = list(MNAR = coef.list1, MAR = coef.list2, MCAR = coef.list3)

###### option2: prevalence
  # get prevalence"
  prev = get_prevalence(mylist)
  # MNAR:
  coef.list1 = vector("list", length = length(var_name))
  for(j in 1:length(var_name)){
    dat.list = vector("list", length = K)
    for(i in 1:(K-K0)){
      dat = mylist$dat.imp.mice.local.option2$MNAR[[i]]$dat
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
      dat = mylist$dat.imp.mice.local.option2$MAR[[i]]$dat
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
      dat = mylist$dat.imp.mice.local.option2$MCAR[[i]]$dat
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
  mylist$coef_DAC_mice_option2 = list(MNAR = coef.list1, MAR = coef.list2, MCAR = coef.list3)
  return(mylist)
}


imp.mice.DAC.all <- function(mylist){
  K = mylist$K
  K0 = K- length(mylist$site_miss) # of sites with no missingness
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
    dat_new1 = mylist$dat.imp.mice.local.option1$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2 = mylist$dat.imp.mice.local.option1$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3 = mylist$dat.imp.mice.local.option1$MCAR[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # MNAR:
      coef = mylist$coef_DAC_mice_option1$MNAR[[j]]
      dat_new1 = imp.mice.pool.1(dat_original1, dat_new1, coef, var_name[j], family[j])
      # MAR:
      coef = mylist$coef_DAC_mice_option1$MAR[[j]]
      dat_new2 = imp.mice.pool.1(dat_original2, dat_new2, coef, var_name[j], family[j])
      # MCAR:
      coef = mylist$coef_DAC_mice_option1$MCAR[[j]]
      dat_new3 = imp.mice.pool.1(dat_original3, dat_new3, coef, var_name[j], family[j]) 
      
    }
    dat_MNAR[[i]] = dat_new1
    dat_MAR[[i]] = dat_new2
    dat_MCAR[[i]] = dat_new3

  }
  mylist$dat.imp.DAC.mice.option1 = list(MNAR = dat_MNAR, MAR = dat_MAR, MCAR = dat_MCAR)
# ####### option2: prevalence
  dat_MNAR = vector("list", length = K-K0)
  dat_MAR = vector("list", length = K-K0)
  dat_MCAR = vector("list", length = K-K0)
  # get prevalence:
  prev = get_prevalence(mylist)
  # loop through each site:
  for(i in 1:(K-K0)){
    # MNAR:
    dat_original1 = mylist$dat_miss$MNAR[[i]] # original missingness
    dat_new1 = mylist$dat.imp.mice.local.option2$MNAR[[i]]$dat
    # MAR:
    dat_original2 = mylist$dat_miss$MAR[[i]] # original missingness
    dat_new2 = mylist$dat.imp.mice.local.option2$MAR[[i]]$dat
    # MCAR:
    dat_original3 = mylist$dat_miss$MCAR[[i]] # original missingness
    dat_new3 = mylist$dat.imp.mice.local.option2$MCAR[[i]]$dat
    # loop through all missing variables:
    for(j in 1:length(var_name)){
      # MNAR:
      coef = mylist$coef_DAC_mice_option2$MNAR[[j]]
      dat_new1 = imp.mice.pool.2(dat_original1, dat_new1, coef, var_name[j], family[j], prev[j])
      # MAR:
      coef = mylist$coef_DAC_mice_option2$MAR[[j]]
      dat_new2 = imp.mice.pool.2(dat_original2, dat_new2, coef, var_name[j], family[j], prev[j])
      # MCAR:
      coef = mylist$coef_DAC_mice_option2$MCAR[[j]]
      dat_new3= imp.mice.pool.2(dat_original3, dat_new3, coef, var_name[j], family[j], prev[j])

    }
    dat_MNAR[[i]] = dat_new1
    dat_MAR[[i]] = dat_new2
    dat_MCAR[[i]] = dat_new3

  }
  mylist$dat.imp.DAC.mice.option2 = list(MNAR = dat_MNAR, MAR = dat_MAR, MCAR = dat_MCAR)
  return(mylist)
}















