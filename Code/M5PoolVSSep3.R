library(tidyverse)
library(gurobi)
library(nlshrink)
library(MASS)

load("Result/MeanPredAll.Rdata")
load("Result/QtlPred.Rdata")
source("Code/func.R")
source("CodePool/VillePool.r")

# 22 and 23 are identical
n_experts = 50
n_sub = 15
load(paste0("ResultPool/MaxMinIdx",n_sub,".Rdata"))

n_days = 28

cov_type = "EW"

min_idx$L2 = 1:n_sub
min_idx$L3 = 1:n_sub
min_idx$L4 = 1:n_sub
min_idx$L5 = 1:n_sub
min_idx$L6 = 1:n_sub
min_idx$L7 = 1:n_sub

load(paste0("ResultPool/M5_",n_sub,".Rdata"))
rm(list = c("display_mat_L7","display_mat_L8","display_mat_L9"))

# L2 ----------------------------------------------------------------------
if(!exists("display_mat_L2")){
  f_data = pred_L2_all
  true_data = df_true_L2
  
  scale2 = scale2_L2[,2]
  weight_vec = (df_weight %>% filter(Level_id == "Level2"))$weight 
  
  f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
  true_data[,-1] = true_data[,-1]/scale2
  
  f_comb = f_data %>% arrange(state_id, group_id) %>% dplyr::select(-state_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
  
  n_state = length(state_vec)
  
  min_pool_idx = rep(min_idx$L2, n_state) + rep(0:(n_state - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L2, n_state) + rep(0:(n_state - 1) * n_experts, each = n_sub)
  
  # u_pool_L2_min_sample = getPooledU(min_pool_idx, "Sample", n_state, n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L2_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_state, n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L2_min_Linear = getPooledU(min_pool_idx, "Linear", n_state, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L2_min_Ville = getPooledU(min_pool_idx, "Ville", n_state, n_sub,
                                   f_comb, u_comb, true_data)
  
  # u_pool_L2_max_sample = getPooledU(max_pool_idx, "Sample", n_state, n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L2_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_state, n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L2_max_Linear = getPooledU(max_pool_idx, "Linear", n_state, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L2_max_Ville = getPooledU(max_pool_idx, "Ville", n_state, n_sub,
                                   f_comb, u_comb, true_data)
  
  u_sep_L2_min_sample = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_min_NonLin = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_min_Linear = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_max_sample = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_max_NonLin = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_max_Linear = matrix(NA, nrow = n_state, ncol = n_days)
  u_EW_L2_min = matrix(NA, nrow = n_state, ncol = n_days)
  u_EW_L2_max = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_min_MS = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_max_MS = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_min_Rob = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_max_Rob = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_min_Cor = matrix(NA, nrow = n_state, ncol = n_days)
  u_sep_L2_max_Cor = matrix(NA, nrow = n_state, ncol = n_days)
  
  
  for(i in 1:n_state){
    f_sep = f_data %>% filter(state_id == state_vec[i]) %>% dplyr::select(-state_id, -group_id)
    true_sep = true_data %>% filter(state_id == state_vec[i]) %>% dplyr::select(-state_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    # group_idx = 1:n_days + (i-1) * n_days
    
    u_sep_L2_min_sample[i,] = unlist(getSepU(min_idx$L2, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L2_min_NonLin[i,] = unlist(getSepU(min_idx$L2, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L2_min_Linear[i,] = unlist(getSepU(min_idx$L2, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L2_min[i,] = unlist(colMeans(u_sep[min_idx$L2,]))
    u_sep_L2_min_MS[i,] = unlist(getSepU(min_idx$L2, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L2_min_Rob[i,] = unlist(getSepU(min_idx$L2, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L2_min_Cor[i,] = unlist(getSepU(min_idx$L2, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    
    u_sep_L2_max_sample[i,] = unlist(getSepU(max_idx$L2, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L2_max_NonLin[i,] = unlist(getSepU(max_idx$L2, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L2_max_Linear[i,] = unlist(getSepU(max_idx$L2, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L2_max[i,] = unlist(colMeans(u_sep[max_idx$L2,]))
    u_sep_L2_max_MS[i,] = unlist(getSepU(max_idx$L2, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L2_max_Rob[i,] = unlist(getSepU(max_idx$L2, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L2_max_Cor[i,] = unlist(getSepU(max_idx$L2, "Cor", n = n_sub, f_sep, u_sep, true_sep))
  }
  
  display_mat_L2 = matrix(c(#getWRMSSE(u_pool_L2_min_sample, weight_vec),
    getWRMSSE(u_pool_L2_min_Linear, weight_vec),
    getWRMSSE(u_pool_L2_min_Ville, weight_vec),
    #getWRMSSE(u_pool_L2_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L2_min_sample, weight_vec),
    getWRMSSE(u_sep_L2_min_Linear, weight_vec),
    getWRMSSE(u_EW_L2_min, weight_vec),
    getWRMSSE(u_sep_L2_min_MS, weight_vec),
    getWRMSSE(u_sep_L2_min_Rob, weight_vec),
    getWRMSSE(u_sep_L2_min_Cor, weight_vec),
    getWRMSSE(u_sep_L2_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L2_max_sample, weight_vec),
    getWRMSSE(u_pool_L2_max_Linear, weight_vec),
    getWRMSSE(u_pool_L2_max_Ville, weight_vec),
    #getWRMSSE(u_pool_L2_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L2_max_sample, weight_vec),
    getWRMSSE(u_sep_L2_max_Linear, weight_vec),
    getWRMSSE(u_EW_L2_max, weight_vec),
    getWRMSSE(u_sep_L2_max_MS, weight_vec),
    getWRMSSE(u_sep_L2_max_Rob, weight_vec),
    getWRMSSE(u_sep_L2_max_Cor, weight_vec),
    getWRMSSE(u_sep_L2_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L2) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L2) = c("min", "max")
  
  print(round(display_mat_L2, 3))
}



# L4 ----------------------------------------------------------------------
if(!exists("display_mat_L4")){
  f_data = pred_L4_all
  true_data = df_true_L4
  
  scale2 = scale2_L4[,2]
  weight_vec = (df_weight %>% filter(Level_id == "Level4"))$weight 
  
  f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
  true_data[,-1] = true_data[,-1]/scale2
  
  f_comb = f_data %>% arrange(cat_id, group_id) %>% dplyr::select(-cat_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
  
  n_cat = length(cat_vec)
  
  min_pool_idx = rep(min_idx$L4, n_cat) + rep(0:(n_cat - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L4, n_cat) + rep(0:(n_cat - 1) * n_experts, each = n_sub)
  
  # u_pool_L4_min_sample = getPooledU(min_pool_idx, "Sample", n_cat, n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L4_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_cat, n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L4_min_Linear = getPooledU(min_pool_idx, "Linear", n_cat, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L4_min_Ville = getPooledU(min_pool_idx, "Ville", n_cat, n_sub,
                                   f_comb, u_comb, true_data)
  
  # u_pool_L4_max_sample = getPooledU(max_pool_idx, "Sample", n_cat, n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L4_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_cat, n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L4_max_Linear = getPooledU(max_pool_idx, "Linear", n_cat, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L4_max_Ville = getPooledU(max_pool_idx, "Ville", n_cat, n_sub,
                                   f_comb, u_comb, true_data)
  
  
  u_sep_L4_min_sample = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_min_NonLin = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_min_Linear = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_max_sample = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_max_NonLin = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_max_Linear = matrix(NA, nrow = n_cat, ncol = n_days)
  u_EW_L4_min = matrix(NA, nrow = n_cat, ncol = n_days)
  u_EW_L4_max = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_min_MS = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_max_MS = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_min_Rob = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_max_Rob = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_min_Cor = matrix(NA, nrow = n_cat, ncol = n_days)
  u_sep_L4_max_Cor = matrix(NA, nrow = n_cat, ncol = n_days)
  
  
  for(i in 1:n_cat){
    f_sep = f_data %>% filter(cat_id == cat_vec[i]) %>% dplyr::select(-cat_id, -group_id)
    true_sep = true_data %>% filter(cat_id == cat_vec[i]) %>% dplyr::select(-cat_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    # group_idx = 1:n_days + (i-1) * n_days
    
    u_sep_L4_min_sample[i,] = unlist(getSepU(min_idx$L4, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L4_min_NonLin[i,] = unlist(getSepU(min_idx$L4, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L4_min_Linear[i,] = unlist(getSepU(min_idx$L4, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L4_min[i,] = unlist(colMeans(u_sep[min_idx$L4,]))
    u_sep_L4_min_MS[i,] = unlist(getSepU(min_idx$L4, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L4_min_Rob[i,] = unlist(getSepU(min_idx$L4, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L4_min_Cor[i,] = unlist(getSepU(min_idx$L4, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    
    u_sep_L4_max_sample[i,] = unlist(getSepU(max_idx$L4, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L4_max_NonLin[i,] = unlist(getSepU(max_idx$L4, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L4_max_Linear[i,] = unlist(getSepU(max_idx$L4, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L4_max[i,] = unlist(colMeans(u_sep[max_idx$L4,]))
    u_sep_L4_max_MS[i,] = unlist(getSepU(max_idx$L4, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L4_max_Rob[i,] = unlist(getSepU(max_idx$L4, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L4_max_Cor[i,] = unlist(getSepU(max_idx$L4, "Cor", n = n_sub, f_sep, u_sep, true_sep))
  }
  
  # u_pool_L4_max_Ville = u_sep_L4_min_NonLin
  
  display_mat_L4 = matrix(c(#getWRMSSE(u_pool_L4_min_sample, weight_vec),
    getWRMSSE(u_pool_L4_min_Linear, weight_vec),
    getWRMSSE(u_pool_L4_min_Ville, weight_vec),
    # getWRMSSE(u_pool_L4_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L4_min_sample, weight_vec),
    getWRMSSE(u_sep_L4_min_Linear, weight_vec),
    getWRMSSE(u_EW_L4_min, weight_vec),
    getWRMSSE(u_sep_L4_min_MS, weight_vec),
    getWRMSSE(u_sep_L4_min_Rob, weight_vec),
    getWRMSSE(u_sep_L4_min_Cor, weight_vec),
    getWRMSSE(u_sep_L4_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L4_max_sample, weight_vec),
    getWRMSSE(u_pool_L4_max_Linear, weight_vec),
    getWRMSSE(u_pool_L4_max_Ville, weight_vec),
    #getWRMSSE(u_pool_L4_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L4_max_sample, weight_vec),
    getWRMSSE(u_sep_L4_max_Linear, weight_vec),
    getWRMSSE(u_EW_L4_max, weight_vec),
    getWRMSSE(u_sep_L4_max_MS, weight_vec),
    getWRMSSE(u_sep_L4_max_Rob, weight_vec),
    getWRMSSE(u_sep_L4_max_Cor, weight_vec),
    getWRMSSE(u_sep_L4_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L4) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L4) = c("min", "max")
  
  print(round(display_mat_L4, 3))
}

# L3 ----------------------------------------------------------------------
if(!exists("display_mat_L3")){
  f_data = pred_L3_all
  true_data = df_true_L3
  
  scale2 = scale2_L3[,2]
  weight_vec = (df_weight %>% filter(Level_id == "Level3"))$weight 
  
  f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
  true_data[,-1] = true_data[,-1]/scale2
  
  f_comb = f_data %>% arrange(store_id, group_id) %>% dplyr::select(-store_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
  
  n_store = length(store_vec)
  
  min_pool_idx = rep(min_idx$L3, n_store) + rep(0:(n_store - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L3, n_store) + rep(0:(n_store - 1) * n_experts, each = n_sub)
  
  # u_pool_L3_min_sample = getPooledU(min_pool_idx, "Sample", n_store, n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L3_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_store, n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L3_min_Linear = getPooledU(min_pool_idx, "Linear", n_store, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L3_min_Ville = getPooledU(min_pool_idx, "Ville", n_store, n_sub,
                                   f_comb, u_comb, true_data)
  
  # u_pool_L3_max_sample = getPooledU(max_pool_idx, "Sample", n_store, n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L3_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_store, n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L3_max_Linear = getPooledU(max_pool_idx, "Linear", n_store, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L3_max_Ville = getPooledU(max_pool_idx, "Ville", n_store, n_sub,
                                   f_comb, u_comb, true_data)
  
  u_sep_L3_min_sample = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_min_NonLin = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_min_Linear = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_max_sample = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_max_NonLin = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_max_Linear = matrix(NA, nrow = n_store, ncol = n_days)
  u_EW_L3_min = matrix(NA, nrow = n_store, ncol = n_days)
  u_EW_L3_max = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_min_MS = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_max_MS = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_min_Rob = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_max_Rob = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_min_Cor = matrix(NA, nrow = n_store, ncol = n_days)
  u_sep_L3_max_Cor = matrix(NA, nrow = n_store, ncol = n_days)
  
  
  for(i in 1:n_store){
    f_sep = f_data %>% filter(store_id == store_vec[i]) %>% dplyr::select(-store_id, -group_id)
    true_sep = true_data %>% filter(store_id == store_vec[i]) %>% dplyr::select(-store_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    # group_idx = 1:n_days + (i-1) * n_days
    
    u_sep_L3_min_sample[i,] = unlist(getSepU(min_idx$L3, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L3_min_NonLin[i,] = unlist(getSepU(min_idx$L3, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L3_min_Linear[i,] = unlist(getSepU(min_idx$L3, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L3_min[i,] = unlist(colMeans(u_sep[min_idx$L3,]))
    u_sep_L3_min_MS[i,] = unlist(getSepU(min_idx$L3, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L3_min_Rob[i,] = unlist(getSepU(min_idx$L3, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L3_min_Cor[i,] = unlist(getSepU(min_idx$L3, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    
    u_sep_L3_max_sample[i,] = unlist(getSepU(max_idx$L3, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L3_max_NonLin[i,] = unlist(getSepU(max_idx$L3, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L3_max_Linear[i,] = unlist(getSepU(max_idx$L3, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L3_max[i,] = unlist(colMeans(u_sep[max_idx$L3,]))
    u_sep_L3_max_MS[i,] = unlist(getSepU(max_idx$L3, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L3_max_Rob[i,] = unlist(getSepU(max_idx$L3, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L3_max_Cor[i,] = unlist(getSepU(max_idx$L3, "Cor", n = n_sub, f_sep, u_sep, true_sep))
  }
  
  display_mat_L3 = matrix(c(#getWRMSSE(u_pool_L3_min_sample, weight_vec),
    getWRMSSE(u_pool_L3_min_Linear, weight_vec),
    getWRMSSE(u_pool_L3_min_Ville, weight_vec),
    #getWRMSSE(u_pool_L3_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L3_min_sample, weight_vec),
    getWRMSSE(u_sep_L3_min_Linear, weight_vec),
    getWRMSSE(u_EW_L3_min, weight_vec),
    getWRMSSE(u_sep_L3_min_MS, weight_vec),
    getWRMSSE(u_sep_L3_min_Rob, weight_vec),
    getWRMSSE(u_sep_L3_min_Cor, weight_vec),
    getWRMSSE(u_sep_L3_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L3_max_sample, weight_vec),
    getWRMSSE(u_pool_L3_max_Linear, weight_vec),
    getWRMSSE(u_pool_L3_max_Ville, weight_vec),
    #getWRMSSE(u_pool_L3_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L3_max_sample, weight_vec),
    getWRMSSE(u_sep_L3_max_Linear, weight_vec),
    getWRMSSE(u_EW_L3_max, weight_vec),
    getWRMSSE(u_sep_L3_max_MS, weight_vec),
    getWRMSSE(u_sep_L3_max_Rob, weight_vec),
    getWRMSSE(u_sep_L3_max_Cor, weight_vec),
    getWRMSSE(u_sep_L3_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L3) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L3) = c("min", "max")
  
  print(round(display_mat_L3, 3))
}

# L5 ----------------------------------------------------------------------
if(!exists("display_mat_L5")){
  f_data = pred_L5_all
  true_data = df_true_L5
  
  scale2 = scale2_L5[,2]
  weight_vec = (df_weight %>% filter(Level_id == "Level5"))$weight 
  
  f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
  true_data[,-1] = true_data[,-1]/scale2
  
  f_comb = f_data %>% arrange(dept_id, group_id) %>% dplyr::select(-dept_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
  
  n_dept = length(dept_vec)
  
  min_pool_idx = rep(min_idx$L5, n_dept) + rep(0:(n_dept - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L5, n_dept) + rep(0:(n_dept - 1) * n_experts, each = n_sub)
  
  # u_pool_L5_min_sample = getPooledU(min_pool_idx, "Sample", n_dept, n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L5_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_dept, n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L5_min_Linear = getPooledU(min_pool_idx, "Linear", n_dept, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L5_min_Ville = getPooledU(min_pool_idx, "Ville", n_dept, n_sub,
                                   f_comb, u_comb, true_data)
  
  # u_pool_L5_max_sample = getPooledU(max_pool_idx, "Sample", n_dept, n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L5_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_dept, n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L5_max_Linear = getPooledU(max_pool_idx, "Linear", n_dept, n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L5_max_Ville = getPooledU(max_pool_idx, "Ville", n_dept, n_sub,
                                   f_comb, u_comb, true_data)
  
  u_sep_L5_min_sample = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_min_NonLin = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_min_Linear = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_max_sample = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_max_NonLin = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_max_Linear = matrix(NA, nrow = n_dept, ncol = n_days)
  u_EW_L5_min = matrix(NA, nrow = n_dept, ncol = n_days)
  u_EW_L5_max = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_min_MS = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_max_MS = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_min_Rob = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_max_Rob = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_min_Cor = matrix(NA, nrow = n_dept, ncol = n_days)
  u_sep_L5_max_Cor = matrix(NA, nrow = n_dept, ncol = n_days)
  
  
  for(i in 1:n_dept){
    f_sep = f_data %>% filter(dept_id == dept_vec[i]) %>% dplyr::select(-dept_id, -group_id)
    true_sep = true_data %>% filter(dept_id == dept_vec[i]) %>% dplyr::select(-dept_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    # group_idx = 1:n_days + (i-1) * n_days
    
    u_sep_L5_min_sample[i,] = unlist(getSepU(min_idx$L5, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L5_min_NonLin[i,] = unlist(getSepU(min_idx$L5, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L5_min_Linear[i,] = unlist(getSepU(min_idx$L5, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L5_min[i,] = unlist(colMeans(u_sep[min_idx$L5,]))
    u_sep_L5_min_MS[i,] = unlist(getSepU(min_idx$L5, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L5_min_Rob[i,] = unlist(getSepU(min_idx$L5, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L5_min_Cor[i,] = unlist(getSepU(min_idx$L5, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    
    u_sep_L5_max_sample[i,] = unlist(getSepU(max_idx$L5, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    # u_sep_L5_max_NonLin[i,] = unlist(getSepU(max_idx$L5, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L5_max_Linear[i,] = unlist(getSepU(max_idx$L5, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_L5_max[i,] = unlist(colMeans(u_sep[max_idx$L5,]))
    u_sep_L5_max_MS[i,] = unlist(getSepU(max_idx$L5, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L5_max_Rob[i,] = unlist(getSepU(max_idx$L5, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_L5_max_Cor[i,] = unlist(getSepU(max_idx$L5, "Cor", n = n_sub, f_sep, u_sep, true_sep))
  }
  
  display_mat_L5 = matrix(c(#getWRMSSE(u_pool_L5_min_sample, weight_vec),
    getWRMSSE(u_pool_L5_min_Linear, weight_vec),
    getWRMSSE(u_pool_L5_min_Ville, weight_vec),
    #getWRMSSE(u_pool_L5_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L5_min_sample, weight_vec),
    getWRMSSE(u_sep_L5_min_Linear, weight_vec),
    getWRMSSE(u_EW_L5_min, weight_vec),
    getWRMSSE(u_sep_L5_min_MS, weight_vec),
    getWRMSSE(u_sep_L5_min_Rob, weight_vec),
    getWRMSSE(u_sep_L5_min_Cor, weight_vec),
    getWRMSSE(u_sep_L5_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L5_max_sample, weight_vec),
    getWRMSSE(u_pool_L5_max_Linear, weight_vec),
    getWRMSSE(u_pool_L5_max_Ville, weight_vec),
    # getWRMSSE(u_pool_L5_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L5_max_sample, weight_vec),
    getWRMSSE(u_sep_L5_max_Linear, weight_vec),
    getWRMSSE(u_EW_L5_max, weight_vec),
    getWRMSSE(u_sep_L5_max_MS, weight_vec),
    getWRMSSE(u_sep_L5_max_Rob, weight_vec),
    getWRMSSE(u_sep_L5_max_Cor, weight_vec),
    getWRMSSE(u_sep_L5_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L5) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L5) = c("min", "max")
  
  print(round(display_mat_L5, 3))
}

# L6 ----------------------------------------------------------------------
if(!exists("display_mat_L6")){
  f_data = pred_L6_all
  true_data = df_true_L6
  
  scale2 = unlist(scale2_L6[,3])
  weight_vec = (df_weight %>% filter(Level_id == "Level6"))$weight 
  
  f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
  true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
  
  n_state = length(state_vec)
  n_cat = length(cat_vec)
  
  min_pool_idx = rep(min_idx$L6, n_state * n_cat) + 
    rep(0:(n_state * n_cat - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L6, n_state * n_cat) +
    rep(0:(n_state *n_cat - 1) * n_experts, each = n_sub)
  
  f_comb = f_data %>% arrange(state_id, cat_id, group_id) %>% ungroup() %>% 
    dplyr::select(-state_id, - cat_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
  
  # u_pool_L6_min_sample = getPooledU(min_pool_idx, "Sample", n_state * n_cat, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L6_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_state * n_cat, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L6_min_Linear = getPooledU(min_pool_idx, "Linear", n_state * n_cat, n = n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L6_min_Ville = getPooledU(min_pool_idx, "Ville", n_state * n_cat, n = n_sub,
                                    f_comb, u_comb, true_data)
  
  # u_pool_L6_max_sample = getPooledU(max_pool_idx, "Sample", n_state * n_cat, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L6_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_state * n_cat, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L6_max_Linear = getPooledU(max_pool_idx, "Linear", n_state * n_cat, n = n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L6_max_Ville = getPooledU(max_pool_idx, "Ville", n_state * n_cat, n = n_sub,
                                    f_comb, u_comb, true_data)
  
  
  u_sep_L6_min_sample = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_min_NonLin = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_min_Linear = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_max_sample = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_max_NonLin = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_max_Linear = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_EW_L6_min = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_EW_L6_max = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_min_MS = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_max_MS = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_min_Rob = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_max_Rob = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_min_Cor = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  u_sep_L6_max_Cor = matrix(NA, nrow = n_state * n_cat, ncol = n_days)
  
  
  for(i in 1:n_state){
    for(j in 1:n_cat){
      idx = (i-1)*n_cat + j
      f_sep = f_data %>% filter(state_id == state_vec[i], cat_id == cat_vec[j]) %>% 
        ungroup() %>% dplyr::select(-state_id, -cat_id, -group_id)
      true_sep = true_data %>% filter(state_id == state_vec[i], cat_id == cat_vec[j]) %>% 
        ungroup() %>% dplyr::select(-state_id, -cat_id)
      
      u_sep = t(t(f_sep) - as.numeric(true_sep))
      
      # group_idx = 1:n_days + (i-1) * n_days
      
      u_sep_L6_min_sample[idx,] = unlist(getSepU(min_idx$L6, "Sample", n = n_sub, f_sep, u_sep, true_sep))
      # u_sep_L6_min_NonLin[idx,] = unlist(getSepU(min_idx$L6, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L6_min_Linear[idx,] = unlist(getSepU(min_idx$L6, "Linear", n = n_sub, f_sep, u_sep, true_sep))
      u_EW_L6_min[idx,] = unlist(colMeans(u_sep[min_idx$L6,]))
      u_sep_L6_min_MS[idx,] = unlist(getSepU(min_idx$L6, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L6_min_Rob[idx,] = unlist(getSepU(min_idx$L6, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L6_min_Cor[idx,] = unlist(getSepU(min_idx$L6, "Cor", n = n_sub, f_sep, u_sep, true_sep))
      
      
      u_sep_L6_max_sample[idx,] = unlist(getSepU(max_idx$L6, "Sample", f_sep, n = n_sub, u_sep, true_sep))
      # u_sep_L6_max_NonLin[idx,] = unlist(getSepU(max_idx$L6, "NonLin", f_sep, n = n_sub, u_sep, true_sep))
      u_sep_L6_max_Linear[idx,] = unlist(getSepU(max_idx$L6, "Linear", f_sep, n = n_sub, u_sep, true_sep))
      u_EW_L6_max[idx,] = unlist(colMeans(u_sep[max_idx$L6,]))
      u_sep_L6_max_MS[idx,] = unlist(getSepU(max_idx$L6, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L6_max_Rob[idx,] = unlist(getSepU(max_idx$L6, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L6_max_Cor[idx,] = unlist(getSepU(max_idx$L6, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    }
  }
  
  # u_pool_L6_max_Ville = u_sep_L6_max_NonLin
  
  display_mat_L6 = matrix(c(#getWRMSSE(u_pool_L6_min_sample, weight_vec),
    getWRMSSE(u_pool_L6_min_Linear, weight_vec),
    getWRMSSE(u_pool_L6_min_Ville, weight_vec),
    #getWRMSSE(u_pool_L6_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L6_min_sample, weight_vec),
    getWRMSSE(u_sep_L6_min_Linear, weight_vec),
    getWRMSSE(u_EW_L6_min, weight_vec),
    getWRMSSE(u_sep_L6_min_MS, weight_vec),
    getWRMSSE(u_sep_L6_min_Rob, weight_vec),
    getWRMSSE(u_sep_L6_min_Cor, weight_vec),
    getWRMSSE(u_sep_L6_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L6_max_sample, weight_vec),
    getWRMSSE(u_pool_L6_max_Linear, weight_vec),
    getWRMSSE(u_pool_L6_max_Ville, weight_vec),
    #getWRMSSE(u_pool_L6_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L6_max_sample, weight_vec),
    getWRMSSE(u_sep_L6_max_Linear, weight_vec),
    getWRMSSE(u_EW_L6_max, weight_vec),
    getWRMSSE(u_sep_L6_max_MS, weight_vec),
    getWRMSSE(u_sep_L6_max_Rob, weight_vec),
    getWRMSSE(u_sep_L6_max_Cor, weight_vec),
    getWRMSSE(u_sep_L6_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L6) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L6) = c("min", "max")
  
  print(round(display_mat_L2, 3))
  print(round(display_mat_L3, 3))
  print(round(display_mat_L4, 3))
  print(round(display_mat_L5, 3))
  print(round(display_mat_L6, 3))
}


# L7 ----------------------------------------------------------------------
if(!exists("display_mat_L7")){
  f_data = pred_L7_all
  true_data = df_true_L7
  
  scale2 = unlist(scale2_L7[,3])
  weight_vec = (df_weight %>% filter(Level_id == "Level7"))$weight
  
  f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
  true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
  
  n_state = length(state_vec)
  n_dept = length(dept_vec)
  
  min_pool_idx = rep(min_idx$L7, n_state * n_dept) +
    rep(0:(n_state * n_dept - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L7, n_state * n_dept) +
    rep(0:(n_state *n_dept - 1) * n_experts, each = n_sub)
  
  f_comb = f_data %>% arrange(state_id, dept_id, group_id) %>% ungroup() %>%
    dplyr::select(-state_id, - dept_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
  
  # u_pool_L7_min_sample = getPooledU(min_pool_idx, "Sample", n_state * n_dept, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L7_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_state * n_dept, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L7_min_Linear = getPooledU(min_pool_idx, "Linear", n_state * n_dept, n = n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L7_min_Ville = getPooledU(min_pool_idx, "Ville", n_state * n_dept, n = n_sub,
                                   f_comb, u_comb, true_data)
  
  # u_pool_L7_max_sample = getPooledU(max_pool_idx, "Sample", n_state * n_dept, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L7_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_state * n_dept, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L7_max_Linear = getPooledU(max_pool_idx, "Linear", n_state * n_dept, n = n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L7_max_Ville = getPooledU(max_pool_idx, "Ville", n_state * n_dept, n = n_sub,
                                   f_comb, u_comb, true_data)
  
  u_sep_L7_min_sample = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_min_NonLin = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_min_Linear = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_max_sample = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_max_NonLin = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_max_Linear = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_EW_L7_min = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_EW_L7_max = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_min_MS = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_max_MS = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_min_Rob = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_max_Rob = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_min_Cor = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  u_sep_L7_max_Cor = matrix(NA, nrow = n_state * n_dept, ncol = n_days)
  
  
  for(i in 1:n_state){
    for(j in 1:n_dept){
      idx = (i-1)*n_dept + j
      f_sep = f_data %>% filter(state_id == state_vec[i], dept_id == dept_vec[j]) %>%
        ungroup() %>% dplyr::select(-state_id, -dept_id, -group_id)
      true_sep = true_data %>% filter(state_id == state_vec[i], dept_id == dept_vec[j]) %>%
        ungroup() %>% dplyr::select(-state_id, -dept_id)
  
      u_sep = t(t(f_sep) - as.numeric(true_sep))
  
      # group_idx = 1:n_days + (i-1) * n_days
  
      u_sep_L7_min_sample[idx,] = unlist(getSepU(min_idx$L7, "Sample", n = n_sub, f_sep, u_sep, true_sep))
      # u_sep_L7_min_NonLin[idx,] = unlist(getSepU(min_idx$L7, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L7_min_Linear[idx,] = unlist(getSepU(min_idx$L7, "Linear", n = n_sub, f_sep, u_sep, true_sep))
      u_EW_L7_min[idx,] = unlist(colMeans(u_sep[min_idx$L7,]))
      u_sep_L7_min_MS[idx,] = unlist(getSepU(min_idx$L7, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L7_min_Rob[idx,] = unlist(getSepU(min_idx$L7, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L7_min_Cor[idx,] = unlist(getSepU(min_idx$L7, "Cor", n = n_sub, f_sep, u_sep, true_sep))
  
  
      u_sep_L7_max_sample[idx,] = unlist(getSepU(max_idx$L7, "Sample", f_sep, n = n_sub, u_sep, true_sep))
      # u_sep_L7_max_NonLin[idx,] = unlist(getSepU(max_idx$L7, "NonLin", f_sep, n = n_sub, u_sep, true_sep))
      u_sep_L7_max_Linear[idx,] = unlist(getSepU(max_idx$L7, "Linear", f_sep, n = n_sub, u_sep, true_sep))
      u_EW_L7_max[idx,] = unlist(colMeans(u_sep[max_idx$L7,]))
      u_sep_L7_max_MS[idx,] = unlist(getSepU(max_idx$L7, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L7_max_Rob[idx,] = unlist(getSepU(max_idx$L7, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L7_max_Cor[idx,] = unlist(getSepU(max_idx$L7, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    }
  }
  
  # u_pool_L7_max_Ville = u_sep_L7_max_NonLin
  
  display_mat_L7 = matrix(c(#getWRMSSE(u_pool_L7_min_sample, weight_vec),
    getWRMSSE(u_pool_L7_min_Linear, weight_vec),
    getWRMSSE(u_pool_L7_min_Ville, weight_vec),
    #getWRMSSE(u_pool_L7_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L7_min_sample, weight_vec),
    getWRMSSE(u_sep_L7_min_Linear, weight_vec),
    getWRMSSE(u_EW_L7_min, weight_vec),
    getWRMSSE(u_sep_L7_min_MS, weight_vec),
    getWRMSSE(u_sep_L7_min_Rob, weight_vec),
    getWRMSSE(u_sep_L7_min_Cor, weight_vec),
    getWRMSSE(u_sep_L7_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L7_max_sample, weight_vec),
    getWRMSSE(u_pool_L7_max_Linear, weight_vec),
    getWRMSSE(u_pool_L7_max_Ville, weight_vec),
    #getWRMSSE(u_pool_L7_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L7_max_sample, weight_vec),
    getWRMSSE(u_sep_L7_max_Linear, weight_vec),
    getWRMSSE(u_EW_L7_max, weight_vec),
    getWRMSSE(u_sep_L7_max_MS, weight_vec),
    getWRMSSE(u_sep_L7_max_Rob, weight_vec),
    getWRMSSE(u_sep_L7_max_Cor, weight_vec),
    getWRMSSE(u_sep_L7_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L7) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L7) = c("min", "max")
  
  print(round(display_mat_L7, 3))
}

# L8 ----------------------------------------------------------------------
if(!exists("display_mat_L8")){
  f_data = pred_L8_all
  true_data = df_true_L8
  
  scale2 = unlist(scale2_L8[,3])
  weight_vec = (df_weight %>% filter(Level_id == "Level8"))$weight
  
  f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
  true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
  
  n_store = length(store_vec)
  n_cat = length(cat_vec)
  
  min_pool_idx = rep(min_idx$L8, n_store * n_cat) +
    rep(0:(n_store * n_cat - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L8, n_store * n_cat) +
    rep(0:(n_store *n_cat - 1) * n_experts, each = n_sub)
  
  f_comb = f_data %>% arrange(store_id, cat_id, group_id) %>% ungroup() %>%
    dplyr::select(-store_id, - cat_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
  
  # u_pool_L8_min_sample = getPooledU(min_pool_idx, "Sample", n_store * n_cat, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L8_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_store * n_cat, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L8_min_Linear = getPooledU(min_pool_idx, "Linear", n_store * n_cat, n = n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L8_min_Ville = getPooledU(min_pool_idx, "Ville", n_store * n_cat, n = n_sub,
                                   f_comb, u_comb, true_data)
  
  # u_pool_L8_max_sample = getPooledU(max_pool_idx, "Sample", n_store * n_cat, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L8_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_store * n_cat, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L8_max_Linear = getPooledU(max_pool_idx, "Linear", n_store * n_cat, n = n_sub,
                                    f_comb, u_comb, true_data)
  # u_pool_L8_max_Ville = getPooledU(max_pool_idx, "Ville", n_store * n_cat, n = n_sub,
  #                                  f_comb, u_comb, true_data)
  u_pool_L8_max_Ville = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  
  u_sep_L8_min_sample = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_min_NonLin = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_min_Linear = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_max_sample = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_max_NonLin = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_max_Linear = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_EW_L8_min = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_EW_L8_max = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_min_MS = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_max_MS = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_min_Rob = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_max_Rob = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_min_Cor = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  u_sep_L8_max_Cor = matrix(NA, nrow = n_store * n_cat, ncol = n_days)
  
  
  for(i in 1:n_store){
    for(j in 1:n_cat){
      idx = (i-1)*n_cat + j
      f_sep = f_data %>% filter(store_id == store_vec[i], cat_id == cat_vec[j]) %>%
        ungroup() %>% dplyr::select(-store_id, -cat_id, -group_id)
      true_sep = true_data %>% filter(store_id == store_vec[i], cat_id == cat_vec[j]) %>%
        ungroup() %>% dplyr::select(-store_id, -cat_id)
      
      u_sep = t(t(f_sep) - as.numeric(true_sep))
      
      # group_idx = 1:n_days + (i-1) * n_days
      
      u_sep_L8_min_sample[idx,] = unlist(getSepU(min_idx$L8, "Sample", n = n_sub, f_sep, u_sep, true_sep))
      # u_sep_L8_min_NonLin[idx,] = unlist(getSepU(min_idx$L8, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L8_min_Linear[idx,] = unlist(getSepU(min_idx$L8, "Linear", n = n_sub, f_sep, u_sep, true_sep))
      u_EW_L8_min[idx,] = unlist(colMeans(u_sep[min_idx$L8,]))
      u_sep_L8_min_MS[idx,] = unlist(getSepU(min_idx$L8, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L8_min_Rob[idx,] = unlist(getSepU(min_idx$L8, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L8_min_Cor[idx,] = unlist(getSepU(min_idx$L8, "Cor", n = n_sub, f_sep, u_sep, true_sep))
      
      
      u_sep_L8_max_sample[idx,] = unlist(getSepU(max_idx$L8, "Sample", f_sep, n = n_sub, u_sep, true_sep))
      # u_sep_L8_max_NonLin[idx,] = unlist(getSepU(max_idx$L8, "NonLin", f_sep, n = n_sub, u_sep, true_sep))
      u_sep_L8_max_Linear[idx,] = unlist(getSepU(max_idx$L8, "Linear", f_sep, n = n_sub, u_sep, true_sep))
      u_EW_L8_max[idx,] = unlist(colMeans(u_sep[max_idx$L8,]))
      u_sep_L8_max_MS[idx,] = unlist(getSepU(max_idx$L8, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L8_max_Rob[idx,] = unlist(getSepU(max_idx$L8, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L8_max_Cor[idx,] = unlist(getSepU(max_idx$L8, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    }
  }
  
  # u_pool_L8_max_Ville = u_sep_L8_max_NonLin
  
  display_mat_L8 = matrix(c(#getWRMSSE(u_pool_L8_min_sample, weight_vec),
    getWRMSSE(u_pool_L8_min_Linear, weight_vec),
    getWRMSSE(u_pool_L8_min_Ville, weight_vec),
    #getWRMSSE(u_pool_L8_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L8_min_sample, weight_vec),
    getWRMSSE(u_sep_L8_min_Linear, weight_vec),
    getWRMSSE(u_EW_L8_min, weight_vec),
    getWRMSSE(u_sep_L8_min_MS, weight_vec),
    getWRMSSE(u_sep_L8_min_Rob, weight_vec),
    getWRMSSE(u_sep_L8_min_Cor, weight_vec),
    getWRMSSE(u_sep_L8_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L8_max_sample, weight_vec),
    getWRMSSE(u_pool_L8_max_Linear, weight_vec),
    getWRMSSE(u_pool_L8_max_Ville, weight_vec),
    #getWRMSSE(u_pool_L8_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L8_max_sample, weight_vec),
    getWRMSSE(u_sep_L8_max_Linear, weight_vec),
    getWRMSSE(u_EW_L8_max, weight_vec),
    getWRMSSE(u_sep_L8_max_MS, weight_vec),
    getWRMSSE(u_sep_L8_max_Rob, weight_vec),
    getWRMSSE(u_sep_L8_max_Cor, weight_vec),
    getWRMSSE(u_sep_L8_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L8) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L8) = c("min", "max")
  
  print(round(display_mat_L8, 3))
}

# L9 ----------------------------------------------------------------------
if(!exists("display_mat_L9")){
  f_data = pred_L9_all
  true_data = df_true_L9
  
  scale2 = unlist(scale2_L9[,3])
  weight_vec = (df_weight %>% filter(Level_id == "Level9"))$weight
  
  f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
  true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
  
  n_store = length(store_vec)
  n_dept = length(dept_vec)
  
  min_pool_idx = rep(min_idx$L9, n_store * n_dept) +
    rep(0:(n_store * n_dept - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx$L9, n_store * n_dept) +
    rep(0:(n_store *n_dept - 1) * n_experts, each = n_sub)
  
  f_comb = f_data %>% arrange(store_id, dept_id, group_id) %>% ungroup() %>%
    dplyr::select(-store_id, - dept_id, -group_id)
  u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
  
  # u_pool_L9_min_sample = getPooledU(min_pool_idx, "Sample", n_store * n_dept, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  # u_pool_L9_min_NonLin = getPooledU(min_pool_idx, "NonLin", n_store * n_dept, n = n_sub,
  #                                      f_comb, u_comb, true_data)
  u_pool_L9_min_Linear = getPooledU(min_pool_idx, "Linear", n_store * n_dept, n = n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L9_min_Ville = getPooledU(min_pool_idx, "Ville", n_store * n_dept, n = n_sub,
                                   f_comb, u_comb, true_data)
  
  # u_pool_L9_max_sample = getPooledU(max_pool_idx, "Sample", n_store * n_dept, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  # u_pool_L9_max_NonLin = getPooledU(max_pool_idx, "NonLin", n_store * n_dept, n = n_sub,
  #                                         f_comb, u_comb, true_data)
  u_pool_L9_max_Linear = getPooledU(max_pool_idx, "Linear", n_store * n_dept, n = n_sub,
                                    f_comb, u_comb, true_data)
  u_pool_L9_max_Ville = getPooledU(max_pool_idx, "Ville", n_store * n_dept, n = n_sub,
                                   f_comb, u_comb, true_data)
  
  u_sep_L9_min_sample = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_min_NonLin = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_min_Linear = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_max_sample = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_max_NonLin = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_max_Linear = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_EW_L9_min = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_EW_L9_max = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_min_MS = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_max_MS = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_min_Rob = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_max_Rob = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_min_Cor = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  u_sep_L9_max_Cor = matrix(NA, nrow = n_store * n_dept, ncol = n_days)
  
  
  for(i in 1:n_store){
    for(j in 1:n_dept){
      idx = (i-1)*n_dept + j
      f_sep = f_data %>% filter(store_id == store_vec[i], dept_id == dept_vec[j]) %>%
        ungroup() %>% dplyr::select(-store_id, -dept_id, -group_id)
      true_sep = true_data %>% filter(store_id == store_vec[i], dept_id == dept_vec[j]) %>%
        ungroup() %>% dplyr::select(-store_id, -dept_id)
      
      u_sep = t(t(f_sep) - as.numeric(true_sep))
      
      # group_idx = 1:n_days + (i-1) * n_days
      
      u_sep_L9_min_sample[idx,] = unlist(getSepU(min_idx$L9, "Sample", n = n_sub, f_sep, u_sep, true_sep))
      # u_sep_L9_min_NonLin[idx,] = unlist(getSepU(min_idx$L9, "NonLin", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L9_min_Linear[idx,] = unlist(getSepU(min_idx$L9, "Linear", n = n_sub, f_sep, u_sep, true_sep))
      u_EW_L9_min[idx,] = unlist(colMeans(u_sep[min_idx$L9,]))
      u_sep_L9_min_MS[idx,] = unlist(getSepU(min_idx$L9, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L9_min_Rob[idx,] = unlist(getSepU(min_idx$L9, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L9_min_Cor[idx,] = unlist(getSepU(min_idx$L9, "Cor", n = n_sub, f_sep, u_sep, true_sep))
      
      
      u_sep_L9_max_sample[idx,] = unlist(getSepU(max_idx$L9, "Sample", f_sep, n = n_sub, u_sep, true_sep))
      # u_sep_L9_max_NonLin[idx,] = unlist(getSepU(max_idx$L9, "NonLin", f_sep, n = n_sub, u_sep, true_sep))
      u_sep_L9_max_Linear[idx,] = unlist(getSepU(max_idx$L9, "Linear", f_sep, n = n_sub, u_sep, true_sep))
      u_EW_L9_max[idx,] = unlist(colMeans(u_sep[max_idx$L9,]))
      u_sep_L9_max_MS[idx,] = unlist(getSepU(max_idx$L9, "MS", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L9_max_Rob[idx,] = unlist(getSepU(max_idx$L9, "Rob", n = n_sub, f_sep, u_sep, true_sep))
      u_sep_L9_max_Cor[idx,] = unlist(getSepU(max_idx$L9, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    }
  }
  
  # u_pool_L9_max_Ville = u_sep_L9_max_NonLin

  display_mat_L9 = matrix(c(#getWRMSSE(u_pool_L9_min_sample, weight_vec),
    getWRMSSE(u_pool_L9_min_Linear, weight_vec),
    getWRMSSE(u_pool_L9_min_Ville, weight_vec),
    #getWRMSSE(u_pool_L9_min_NonLin, weight_vec),
    getWRMSSE(u_sep_L9_min_sample, weight_vec),
    getWRMSSE(u_sep_L9_min_Linear, weight_vec),
    getWRMSSE(u_EW_L9_min, weight_vec),
    getWRMSSE(u_sep_L9_min_MS, weight_vec),
    getWRMSSE(u_sep_L9_min_Rob, weight_vec),
    getWRMSSE(u_sep_L9_min_Cor, weight_vec),
    getWRMSSE(u_sep_L9_min_NonLin, weight_vec),
    #getWRMSSE(u_pool_L9_max_sample, weight_vec),
    getWRMSSE(u_pool_L9_max_Linear, weight_vec),
    getWRMSSE(u_pool_L9_max_Ville, weight_vec),
    #getWRMSSE(u_pool_L9_max_NonLin, weight_vec),
    getWRMSSE(u_sep_L9_max_sample, weight_vec),
    getWRMSSE(u_sep_L9_max_Linear, weight_vec),
    getWRMSSE(u_EW_L9_max, weight_vec),
    getWRMSSE(u_sep_L9_max_MS, weight_vec),
    getWRMSSE(u_sep_L9_max_Rob, weight_vec),
    getWRMSSE(u_sep_L9_max_Cor, weight_vec),
    getWRMSSE(u_sep_L9_max_NonLin, weight_vec)
  ),
  byrow =TRUE, nrow = 2)
  colnames(display_mat_L9) = c("Pool Linear", "Pool Ville","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor", "Sep NonLin")
  rownames(display_mat_L9) = c("min", "max")
  
  print(round(display_mat_L9, 3))
}

save(list = paste0("display_mat_L", 2:9), file = paste0("ResultPool/M5_",n_sub,".Rdata"))

library(xtable)
n_sub = 2
load(paste0("ResultPool/M5_",n_sub,".Rdata"))

display_mat = rbind(display_mat_L2, display_mat_L3, display_mat_L4,
                    display_mat_L5, display_mat_L6, display_mat_L7,
                    display_mat_L8, display_mat_L9)[c(1,3,5,7,9,11,13,15) + 1,c(1,5,3,4,8,6,7)]

# display_mat = rbind(display_mat_L2, display_mat_L3, display_mat_L4,
#                     display_mat_L5, display_mat_L6, display_mat_L7[,-9],
#                     display_mat_L8[,-9], display_mat_L9[,-9])[c(2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15),c(1,2,3,4,8,6,7,5)]


xtable(display_mat, digits = 3)
apply(display_mat, 1, which.min)
