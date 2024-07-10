get_prevalence <- function(mylist){
  K = mylist$K
  S = c(1:K)
  S0 = S[-which(mylist$site_miss %in% S)]
  col_miss = mylist$col_miss
  family = mylist$family
  vec = c()

  for(i in 1:length(col_miss)){
    temp = 0
    for(k in S0){
      # access k-th site with no missingness:
      dat = mylist$dat_original[[k]]$dat
      if(family[i] == "binomial"){
        var = dat[[col_miss[i]]]
        temp = temp + sum(var==1)/length(var)
      }
    }
    vec[i] = ifelse(family[i] == "binomial", temp/length(S0), NA)
    
  }
  return(vec)
}

get_prevalence_realdata <- function(mylist, col_miss, family, site_miss){
  K = length(mylist)
  S = c(1:K)
  S0 = S[-which(site_miss %in% S)]
  vec = c()
  
  for(i in 1:length(col_miss)){
    temp = 0
    for(k in S0){
      dat = mylist[[k]]
      if(family[i] == "binomial"){
        var = dat[[col_miss[i]]]
        temp = temp + sum(var==1)/length(var)
      }
    }
    vec[i] = ifelse(family[i] == "binomial", temp/length(S0), NA)
    
  }
  return(vec)
}

get.missingrate.all <- function(mylist){
  K = mylist$K
  K0 = length(mylist$site_miss)
  col_miss = mylist$col_miss
  J = length(col_miss)
  res = matrix(nrow = 3*K0*J, ncol = 4)
  row = 0
  missing_type = "MNAR"
  for(i in 1:K0){
    for(j in 1:J){
      row = row + 1
      rate_temp = get.missingrate(mylist$dat_miss$MNAR[[i]][[col_miss[j]]])
      res[row,] = c(col_miss[j], rate_temp, paste0("Site",i), missing_type)
    }
  }
  missing_type = "MAR"
  for(i in 1:K0){
    for(j in 1:J){
      row = row + 1
      rate_temp = get.missingrate(mylist$dat_miss$MAR[[i]][[col_miss[j]]])
      res[row,] = c(col_miss[j], rate_temp, paste0("Site",i), missing_type)
    }
  }
  missing_type = "MCAR"
  for(i in 1:K0){
    for(j in 1:J){
      row = row + 1
      rate_temp = get.missingrate(mylist$dat_miss$MCAR[[i]][[col_miss[j]]])
      res[row,] = c(col_miss[j], rate_temp, paste0("Site",i), missing_type)
    }
  }
  res = as.data.frame(res)
  colnames(res) = c("Var", "Missing rate", "Site name", "Missing type")
  return(res)
}

# get coef of non-missingness site for all variables
get.coef.nmiss.all <- function(mylist){
  K = length(mylist$site_miss)
  var_miss = mylist$col_miss
  family = mylist$family
  K0 = mylist$K - K
  J = length(var_miss)
  dat = mylist$dat_original
  coef_nmiss = vector("list", length = K0)
  for(i in 1:K0){
    coef = vector("list", length = J)
    dat_temp = dat[[K+i]]$dat
    for(j in 1:J){
      coef_temp = get.coef.nmiss.single(var_miss[j], dat_temp, family[j])
      coef[[j]] = coef_temp
    }
    names(coef) = var_miss
    coef_nmiss[[i]] = coef
  }
  names(coef_nmiss) = paste0("S", seq(from=K+1, to=K+K0, by=1))
  mylist$coef_non_missing = coef_nmiss
  return(mylist)
}

# get coef of non-missingness site for all variables
get.coef.nmiss.all.new <- function(mylist){
  K = mylist$K
  var_miss = mylist$col_miss
  family = mylist$family
  J = length(var_miss)
  dat = mylist$dat_original
  coef_nmiss = vector("list", length = K)
  for(i in 1:K){
    coef = vector("list", length = J)
    dat_temp = dat[[i]]$dat
    for(j in 1:J){
      coef_temp = get.coef.nmiss.single(var_miss[j], dat_temp, family[j])
      coef[[j]] = coef_temp
    }
    names(coef) = var_miss
    coef_nmiss[[i]] = coef
  }
  names(coef_nmiss) = paste0("S", c(1:K))
  mylist$coef_non_missing = coef_nmiss
  return(mylist)
}

get.coef.nmiss.single <- function(var_name, dat, family){
  dat2 = dat[, -1] #drop response y
  dat3 = subset(dat2, select = -which(names(dat2) == var_name)) %>% as.matrix()
  y =  subset(dat2, select = which(names(dat2) == var_name)) %>% as.matrix()
  res = lasso.fun(x = dat3, y = as.vector(y), family = family)
  coef = res$coef
  return(coef)
}

# check if input of DAC function has NA values
check_NA <- function(list){
  N = length(list)
  for(i in 1:N){
    if(anyNA(list[[i]])){
      print("The input dat.list for DOCS function has NA!")
      print(i)
      return()
    }
  }
  print("No NA detected in the input of DOCS function.")
}

# set a function to create directory for a given directory name
create_dir_mac <- function(parent_dir, sub = c("coef", "dist", "diff", "roc", "summary")){
  if(!file.exists(parent_dir)){dir.create(parent_dir)}
  for (i in 1:length(sub)){
    # Create the full path for the subdirectory
    sub_dir_path <- file.path(parent_dir, sub[i])
    
    # Check if the subdirectory already exists
    if (!file.exists(sub_dir_path)) {
      # Create the subdirectory
      dir.create(sub_dir_path)
      cat("Subdirectory created:", sub_dir_path, "\n")
    } else {
      cat("Subdirectory already exists:", sub_dir_path, "\n")
    }
  }
}

create_subdirectory <- function(parent_dir, sub_dir_name) {
  # Create the full path for the subdirectory
  sub_dir_path <- file.path(parent_dir, sub_dir_name)
  
  # Check if the subdirectory already exists
  if (!file.exists(sub_dir_path)) {
    # Create the subdirectory
    dir.create(sub_dir_path)
    cat("Subdirectory created:", sub_dir_path, "\n")
  } else {
    cat("Subdirectory already exists:", sub_dir_path, "\n")
  }
}

normal_based_CI <- function(est, se, level = 0.05){
  upper <- est + qnorm(1-level/2)*se
  lower <- est - qnorm(1-level/2)*se
  return(cbind(lower, upper))
}

coverage_probability <- function(lower, upper, truth){
  mean(I(upper >= truth)*I(lower <= truth))
}

get_coverage <- function(dat, truth, p, level) {
  N <- nrow(dat)
  se <- lapply(dat, sd)
  coverage <- sapply(1:(p+1), function(j) {
    coverage_sum <- sum(sapply(1:N, function(i) {
      cis <- normal_based_CI(dat[i, j], se[[j]], level)
      coverage_probability(cis[1], cis[2], truth[j])
    }))
    coverage_sum / N
  })
  names(coverage) <- colnames(dat)
  return(coverage)
}

get_bias_prediction<- function(dat, truth, p) {
  N <- nrow(dat)
  bias <- sapply(1:(p+1), function(j) {
    column <- dat[, j]
    abs_diff <- abs(column - truth[j])
    sum(abs_diff) / N
  })
  names(bias) <- colnames(dat)
  return(bias)
}





