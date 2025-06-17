library(mvtnorm)
library(nlshrink)
library(tidyverse)
source("Code/func.R")

set.seed(20240422)

n_experts = 15

n_prod_vec = c(3, 7, 9, 10, 21, 30, 70, 700, 7000) 
length_prod_vec = length(n_prod_vec)

sigma2_x = 1
sigma2_vec = c(rep(1, 7) , rep(5, 8)) #Table 9
# sigma2_vec = c(rep(1, 7) , rep(60, 8)) #Table 10
sigma2_eps = 30

#It takes a long time to run M = 100. One could reduce the time by setting M = 10
M = 100

n_obs = 47
#It takes a long time to run n_test = 3000. One could reduce the time by setting n_test = 1000
n_test = 3000

RMSE_array = array(NA, dim = c(M, 3, length_prod_vec))

time_vec = rep(NA, length_prod_vec)

for(i in 1:length_prod_vec){
  n_products = n_prod_vec[i]
  
  cov_ui = sigma2_x + diag(sigma2_vec + sigma2_eps)
  
  data_test = simData(n_test, sigma2_x, sigma2_vec, sigma2_eps, n_products, n_experts)
  
  start_time = Sys.time()
  
  RMSE_array[,1,i] = 
    map_dbl(1:M, ~ get1SimPoolRMSE(data_test, n_obs, sigma2_x, sigma2_vec, sigma2_eps, n_products, n_experts))
  
  end_time = Sys.time()
  
  time_vec[i] = as.numeric(difftime(end_time, start_time, units = "secs"))
  
  print(paste(n_products, "Avg Time", round(time_vec[i]/M,3), "s"))
  
  RMSE_array[,2,i] = sqrt(mean(cov_ui))
  RMSE_array[,3,i] = sqrt(1/sum(solve(cov_ui)))
}


print(t(apply(RMSE_array, 3, colMeans)), digits = 4)

print(apply(RMSE_array, 3, function(x) mean(x[,1] < x[,3])), digits = 3)


