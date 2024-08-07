library(mvtnorm)
library(nlshrink)
library(purrr)
source("Code/func.R")

set.seed(20240422)

n_experts = 15

n_prod_vec = c(3, 7, 9, 10, 21, 30, 70)
length_prod_vec = length(n_prod_vec)

sigma2_x = 1
# sigma2_vec = c(rep(1, 7) , rep(5, 8)) #Table 9
sigma2_vec = c(rep(1, 7) , rep(60, 8)) #Table 10 
sigma2_eps = 30

#It takes several minutes to run M = 100. One could reduce the time by setting M = 10
M = 1

n_obs = 47

var_array = array(NA, dim = c(M, 3, length_prod_vec))


for(i in 1:length_prod_vec){
  n_products = n_prod_vec[i]
  
  cov_ui = sigma2_x + diag(sigma2_vec + sigma2_eps)
  cov_u = diag(n_products) %x% matrix(sigma2_x, n_experts, n_experts) +
    matrix(1, n_products, n_products) %x% diag(sigma2_vec) +
    sigma2_eps * diag(n_products) %x% diag(n_experts)
  
  var_array[,1,i] = map_dbl(1:M, ~ get1SimPoolVar(cov_u, n_obs, n_products, n_experts))
  
  var_array[,2,i] = mean(cov_ui)
  var_array[,3,i] = 1/sum(solve(cov_ui))
}

print(t(apply(var_array, 3, colMeans)), digits = 4)

print(apply(var_array, 3, function(x) mean(x[,1] < x[,3])), digits = 3)


