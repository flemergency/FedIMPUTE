# impute locally, using mean/median strategy 
imp.local.mean.all <- function(mylist){
  var_list = mylist$col_miss
  family = mylist$family
  site = mylist$site_miss
  # MNAR:
  list1 = vector("list", length = length(site))
  for(i in site){
    dat_temp = mylist$dat_miss$MNAR[[i]]
    dat_temp_new = imp.mean.2(dat_temp, var_list, family)
    list1[[i]] = dat_temp_new
  }
  # MAR:
  list2 = vector("list", length = length(site))
  for(i in site){
    dat_temp = mylist$dat_miss$MAR[[i]]
    dat_temp_new = imp.mean.2(dat_temp, var_list, family)
    list2[[i]] = dat_temp_new
  }
  # MCAR:
  list3 = vector("list", length = length(site))
  for(i in site){
    dat_temp = mylist$dat_miss$MCAR[[i]]
    dat_temp_new = imp.mean.2(dat_temp, var_list, family)
    list3[[i]] = dat_temp_new
  }
  list_temp = list(list1, list2, list3)
  names(list_temp) = c("MNAR", "MAR", "MCAR")
  mylist$dat.imp.local.mean = list_temp
  
  return(mylist)
}

# impute naive-federated, using mean/median from S3 (no-missingness site) for imputation
imp.naive.fed.all <- function(mylist){
  var_list = mylist$col_miss
  family = mylist$family
  site = mylist$site_miss
  N = mylist$N
  NN = N/sum(N)  # normalized weights
  K = mylist$K
  K0 = length(mylist$site_miss)
  res0 = c()
  # first calculate the vec with no-missingness site; 
  # currently only using one site's info
  dat_temp = mylist$dat_original[[K0+1]]$dat
  for(i in 1:length(var_list)){
    var_temp = dat_temp[[var_list[i]]]
    res0_temp = ifelse(family[i]=="binomial",
                       as.numeric(names(sort(table(var_temp), decreasing = TRUE)[1])),
                      mean(var_temp))
    res0 = c(res0, res0_temp)
    }
  # MNAR:
  list1 = vector("list", length = length(site))
  for(i in site){
    dat_temp = mylist$dat_miss$MNAR[[i]]
    dat_temp_new = imp.mean.3(dat_temp, var_list, res0)
    list1[[i]] = dat_temp_new
  }
  # MAR:
  list2 = vector("list", length = length(site))
  for(i in site){
    dat_temp = mylist$dat_miss$MAR[[i]]
    dat_temp_new = imp.mean.3(dat_temp, var_list, res0)
    list2[[i]] = dat_temp_new
  }
  # MCAR:
  list3 = vector("list", length = length(site))
  for(i in site){
    dat_temp = mylist$dat_miss$MCAR[[i]]
    dat_temp_new = imp.mean.3(dat_temp, var_list, res0)
    list3[[i]] = dat_temp_new
  }
  
  list_temp = list(list1, list2, list3)
  names(list_temp) = c("MNAR", "MAR", "MCAR")
  mylist$dat.imp.naive.fed = list_temp
  return(mylist)
}








