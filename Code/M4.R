library(tidyverse)
library(nlshrink)
library(MASS)
source("Code/func.R")

load("CleanedData/M4.Rdata")

n_experts = 17 #total number of experts
n_periods = 48 #number of observations
n_prods = 119

methods = c("EW", "Sample", "Linear", "Cor", "S+EW", "Var", "Rob")
n_methods = length(methods)

display_mat = matrix(NA, nrow = 1, ncol = 1 + length(methods),
                     dimnames = list("All", 
                                     c("P:Linear", "S:EW", "S:Sample", "S:Linear",
                                       "S:Cor", "S:S+EW", "S:Var", "S:Rob")))

true_data = data_test %>% dplyr::select(-V1)
f_comb = f_data %>% dplyr::select(-id, -group_id)
u_comb = u_data %>% dplyr::select(-id, -group_id)

pool_idx = rep(1:n_experts, n_prods) + rep(0:(n_prods - 1) * n_experts, each = n_experts)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_prods, n_experts, f_comb, u_comb, true_data)

u_sep_methods = array(NA, dim = c(n_methods, n_prods, n_periods), dimnames = list(methods))

for(i in 1:n_prods){
  single_prod_idx = 1:n_experts + (i-1) * n_experts
  
  true_sep = true_data[i,]
  f_sep = tibble(f_comb[single_prod_idx,]) 
  u_sep = u_comb[single_prod_idx,]
  
  u_sep_methods[,i,] = getSepUMethods(1:n_experts, methods, n = n_experts, f_sep, u_sep, true_sep)
}

display_mat["All",] = c(getWRMSSE(u_pool_Linear), apply(u_sep_methods, 1, getWRMSSE))

print(display_mat, digits = 3)
## 0.284 1.09 0.453 0.518 0.672 0.431 0.771  0.67