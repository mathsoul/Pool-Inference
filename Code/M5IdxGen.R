library(tidyverse)
library(gurobi)

load("CleanedData/M5.Rdata")
source("Code/func.R")

# 22 and 23 are identical

n_experts = 50
n_top = 17
n_days = 28
pct = 0.1

cov_type = "EW"

max_idx = list()
max_idx_17 = list()
top_idx = list()

for(n_sub in c(2,3,10,15)){
  # L2 ----------------------------------------------------------------------
  if(is.null(max_idx$L2)){
    f_data = pred_L2_all
    true_data = df_true_L2
    
    scale2 = scale2_L2[,2]
    weight_vec = (df_weight %>% filter(Level_id == "Level2"))$weight 
    
    f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
    true_data[,-1] = true_data[,-1]/scale2
    
    n_state = length(state_vec)
    
    f_comb = f_data %>% arrange(state_id, group_id) %>% dplyr::select(-state_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L2 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L2 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L2 = 1:n_sub
    
    stopifnot(length(max_idx$L2) == n_sub & length(top_idx$L2) == n_sub)
  }
  # L3 ----------------------------------------------------------------------
  if(is.null(max_idx$L3)){
    f_data = pred_L3_all
    true_data = df_true_L3
    
    scale2 = scale2_L3[,2]
    weight_vec = (df_weight %>% filter(Level_id == "Level3"))$weight 
    
    f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
    true_data[,-1] = true_data[,-1]/scale2
    
    n_store = length(store_vec)
    
    f_comb = f_data %>% arrange(store_id, group_id) %>% dplyr::select(-store_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L3 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L3 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L3 = 1:n_sub
    
    stopifnot(length(max_idx$L3) == n_sub & length(top_idx$L3) == n_sub)
  }
  
  # L4 ----------------------------------------------------------------------
  if(is.null(max_idx$L4)){
    f_data = pred_L4_all
    true_data = df_true_L4
    
    scale2 = scale2_L4[,2]
    weight_vec = (df_weight %>% filter(Level_id == "Level4"))$weight 
    
    f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
    true_data[,-1] = true_data[,-1]/scale2
    
    f_comb = f_data %>% arrange(cat_id, group_id) %>% dplyr::select(-cat_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L4 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L4 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L4 = 1:n_sub
    
    stopifnot(length(max_idx$L4) == n_sub & length(top_idx$L4) == n_sub)
  }
  
  # L5 ----------------------------------------------------------------------
  if(is.null(max_idx$L5)){
    f_data = pred_L5_all
    true_data = df_true_L5
    
    scale2 = scale2_L5[,2]
    weight_vec = (df_weight %>% filter(Level_id == "Level5"))$weight 
    
    f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
    true_data[,-1] = true_data[,-1]/scale2
    
    f_comb = f_data %>% arrange(dept_id, group_id) %>% dplyr::select(-dept_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-1])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L5 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L5 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L5 = 1:n_sub
    
    stopifnot(length(max_idx$L5) == n_sub & length(top_idx$L5) == n_sub)
  }
  
  # L6 ----------------------------------------------------------------------
  if(is.null(max_idx$L6)){
    f_data = pred_L6_all
    true_data = df_true_L6
    
    scale2 = unlist(scale2_L6[,3])
    weight_vec = (df_weight %>% filter(Level_id == "Level6"))$weight 
    
    f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
    true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
    
    f_comb = f_data %>% arrange(state_id, cat_id, group_id) %>% ungroup() %>% 
      dplyr::select(-state_id, - cat_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L6 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L6 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L6 = 1:n_sub
    
    stopifnot(length(max_idx$L6) == n_sub & length(top_idx$L6) == n_sub)
  }
  
  # L7 ----------------------------------------------------------------------
  if(is.null(max_idx$L7)){
    f_data = pred_L7_all
    true_data = df_true_L7
    
    scale2 = unlist(scale2_L7[,3])
    weight_vec = (df_weight %>% filter(Level_id == "Level7"))$weight 
    
    f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
    true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
    
    f_comb = f_data %>% arrange(state_id, dept_id, group_id) %>% ungroup() %>% 
      dplyr::select(-state_id, - dept_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L7 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L7 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L7 = 1:n_sub
    
    stopifnot(length(max_idx$L7) == n_sub & length(top_idx$L7) == n_sub)
  }
  
  # L8 ----------------------------------------------------------------------
  if(is.null(max_idx$L8)){
    f_data = pred_L8_all
    true_data = df_true_L8
    
    scale2 = unlist(scale2_L8[,3])
    weight_vec = (df_weight %>% filter(Level_id == "Level8"))$weight 
    
    f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
    true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
    
    f_comb = f_data %>% arrange(store_id, cat_id, group_id) %>% ungroup() %>% 
      dplyr::select(-store_id, - cat_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L8 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L8 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L8 = 1:n_sub
    
    stopifnot(length(max_idx$L8) == n_sub & length(top_idx$L8) == n_sub)
  }
  
  # L9 ----------------------------------------------------------------------
  if(is.null(max_idx$L9)){
    f_data = pred_L9_all
    true_data = df_true_L9
    
    scale2 = unlist(scale2_L9[,3])
    weight_vec = (df_weight %>% filter(Level_id == "Level9"))$weight 
    
    f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
    true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2
    
    f_comb = f_data %>% arrange(store_id, dept_id, group_id) %>% ungroup() %>% 
      dplyr::select(-store_id, - dept_id, -group_id)
    u_comb = as.matrix(f_comb - as.matrix(true_data[,-c(1,2)])  %x% rep(1, n_experts))
    
    var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
    
    dist_mat_max = getDistMat(var_comb, weight_vec)
    max_idx$L9 = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max)) == 1)
    max_idx_17$L9 = which(round(OptDivIdx(n_top, n_sub, dist_mat_max[1:n_top,1:n_top])) == 1)
    
    top_idx$L9 = 1:n_sub
    
    stopifnot(length(max_idx$L9) == n_sub & length(top_idx$L9) == n_sub)
  }
  
  idx = top_idx
  save(list = "idx", file = paste0("TeamIndex/TopIdx", n_sub, ".Rdata"))
  
  idx = max_idx
  save(list = "idx", file = paste0("TeamIndex/MaxIdx", n_sub, ".Rdata"))
  
  idx = max_idx_17
  save(list = "idx", file = paste0("TeamIndex/Max17Idx", n_sub, ".Rdata"))
}
