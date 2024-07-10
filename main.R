# demo for real data analysis with FedIMPUTE
source("R/fun.realdata.R")
source("R/fun.impute.single.R")
source("R/fun.impute.mice.R")
source("R/fun.impute.naive.R")
source("R/helpers.R")
source("R/utils.R")

myseed = 123
set.seed(myseed)

# load simulated data:
dat_miss = readRDS("data/mydata.rds")
# specify site with missingness:
site_miss = c(1,2)
# specify var with missingness :
col_miss = c("X1", "X2")
# specify var type:
family = rep("gaussian", 2)


res.all = FedIMPUTE(dat_miss, col_miss, family, site_miss, myseed, maxit_mice = 5, maxit_DAC = 3, ridge = T, ns = F, penalty = T, lam.TL.target=T)
