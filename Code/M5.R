library(tidyverse)
library(nlshrink)
library(MASS)
source("Code/func.R")

load("CleanedData/M5.Rdata")

n_experts = 50 #total number of experts
n_periods = 28 #number of observations

methods = c("EW", "Sample", "Linear", "Cor", "S+EW", "Var", "Rob")
n_methods = length(methods)

team_idx_name = "TeamIndex/TopIdx2.Rdata" #Table 5
# team_idx_name = "TeamIndex/TopIdx15.Rdata" #Table 6
# team_idx_name = "TeamIndex/TopIdx3.Rdata" #Table EC.1
# team_idx_name = "TeamIndex/TopIdx10.Rdata" #Table EC.2
# team_idx_name = "TeamIndex/MaxIdx2.Rdata" #Table EC.5
# team_idx_name = "TeamIndex/MaxIdx15.Rdata" #Table EC.6
# team_idx_name = "TeamIndex/Max17Idx2.Rdata" #Table EC.9
# team_idx_name = "TeamIndex/Max17Idx15.Rdata" #Table EC.10

n_sub = extractNSub(team_idx_name) #number of experts in aggregation
load(team_idx_name)

display_mat = matrix(NA, nrow = 8, ncol = 1 + length(methods),
                     dimnames = list(paste0("L",2:9), 
                                     c("P:Linear", "S:EW", "S:Sample", "S:Linear",
                                       "S:Cor", "S:S+EW", "S:Var", "S:Rob")))


# L2 ----------------------------------------------------------------------
f_data = pred_L2_all
true_data = df_true_L2

scale2 = scale2_L2[,2]
weight_vec = (df_weight %>% filter(Level_id == "Level2"))$weight 

f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
true_data[,-1] = true_data[,-1]/scale2

true_data_pool = true_data %>% dplyr::select(-state_id)
f_pool = f_data %>% arrange(state_id, group_id) %>% dplyr::select(-state_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

n_state = length(state_vec)
  
pool_idx = rep(idx$L2, n_state) + rep(0:(n_state - 1) * n_experts, each = n_sub)
  
u_pool_Linear = getPooledU(pool_idx, "Linear", n_state, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_state, n_periods), dimnames = list(methods))
  
for(i in 1:n_state){
  f_sep = f_data %>% filter(state_id == state_vec[i]) %>% dplyr::select(-state_id, -group_id)
  true_sep = true_data %>% filter(state_id == state_vec[i]) %>% dplyr::select(-state_id)
  
  u_sep = t(t(f_sep) - as.numeric(true_sep))
  
  u_sep_methods[,i,] = getSepUMethods(idx$L2, methods, n = n_sub, f_sep, u_sep, true_sep)
}

display_mat["L2",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))

# L3 ----------------------------------------------------------------------
f_data = pred_L3_all
true_data = df_true_L3

scale2 = scale2_L3[,2]
weight_vec = (df_weight %>% filter(Level_id == "Level3"))$weight 

f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
true_data[,-1] = true_data[,-1]/scale2

true_data_pool = true_data %>% dplyr::select(-store_id)
f_pool = f_data %>% arrange(store_id, group_id) %>% dplyr::select(-store_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

n_store = length(store_vec)

pool_idx = rep(idx$L3, n_store) + rep(0:(n_store - 1) * n_experts, each = n_sub)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_store, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_store, n_periods), dimnames = list(methods))

for(i in 1:n_store){
  f_sep = f_data %>% filter(store_id == store_vec[i]) %>% dplyr::select(-store_id, -group_id)
  true_sep = true_data %>% filter(store_id == store_vec[i]) %>% dplyr::select(-store_id)
  
  u_sep = t(t(f_sep) - as.numeric(true_sep))
  
  u_sep_methods[,i,] = getSepUMethods(idx$L3, methods, n = n_sub, f_sep, u_sep, true_sep)
}

display_mat["L3",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))

# L4 ----------------------------------------------------------------------
f_data = pred_L4_all
true_data = df_true_L4

scale2 = scale2_L4[,2]
weight_vec = (df_weight %>% filter(Level_id == "Level4"))$weight 

f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
true_data[,-1] = true_data[,-1]/scale2

true_data_pool = true_data %>% dplyr::select(-cat_id)
f_pool = f_data %>% arrange(cat_id, group_id) %>% dplyr::select(-cat_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

n_cat = length(cat_vec)

pool_idx = rep(idx$L4, n_cat) + rep(0:(n_cat - 1) * n_experts, each = n_sub)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_cat, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_cat, n_periods), dimnames = list(methods))

for(i in 1:n_cat){
  f_sep = f_data %>% filter(cat_id == cat_vec[i]) %>% dplyr::select(-cat_id, -group_id)
  true_sep = true_data %>% filter(cat_id == cat_vec[i]) %>% dplyr::select(-cat_id)
  
  u_sep = t(t(f_sep) - as.numeric(true_sep))
  
  u_sep_methods[,i,] = getSepUMethods(idx$L4, methods, n = n_sub, f_sep, u_sep, true_sep)
}

display_mat["L4",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))


# L5 ----------------------------------------------------------------------
f_data = pred_L5_all
true_data = df_true_L5

scale2 = scale2_L5[,2]
weight_vec = (df_weight %>% filter(Level_id == "Level5"))$weight 

f_data[,-c(1,30)] = f_data[,-c(1,30)]/scale2
true_data[,-1] = true_data[,-1]/scale2

true_data_pool = true_data %>% dplyr::select(-dept_id)
f_pool = f_data %>% arrange(dept_id, group_id) %>% dplyr::select(-dept_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

n_dept = length(dept_vec)

pool_idx = rep(idx$L5, n_dept) + rep(0:(n_dept - 1) * n_experts, each = n_sub)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_dept, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_dept, n_periods), dimnames = list(methods))

for(i in 1:n_dept){
  f_sep = f_data %>% filter(dept_id == dept_vec[i]) %>% dplyr::select(-dept_id, -group_id)
  true_sep = true_data %>% filter(dept_id == dept_vec[i]) %>% dplyr::select(-dept_id)
  
  u_sep = t(t(f_sep) - as.numeric(true_sep))
  
  u_sep_methods[,i,] = getSepUMethods(idx$L5, methods, n = n_sub, f_sep, u_sep, true_sep)
}

display_mat["L5",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))

# L6 ----------------------------------------------------------------------
f_data = pred_L6_all
true_data = df_true_L6

scale2 = unlist(scale2_L6[,3])
weight_vec = (df_weight %>% filter(Level_id == "Level6"))$weight 

f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2

n_state = length(state_vec)
n_cat = length(cat_vec)

true_data_pool = true_data %>% ungroup() %>%dplyr::select(-state_id, -cat_id)
f_pool = f_data %>% arrange(state_id, cat_id, group_id) %>% ungroup() %>% 
  dplyr::select(-state_id, - cat_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

pool_idx = rep(idx$L6, n_state * n_cat) + 
  rep(0:(n_state * n_cat - 1) * n_experts, each = n_sub)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_state * n_cat, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_state * n_cat, n_periods),
                      dimnames = list(methods))

for(i in 1:n_state){
  for(j in 1:n_cat){
    variable_idx = (i-1)*n_cat + j
    f_sep = f_data %>% filter(state_id == state_vec[i], cat_id == cat_vec[j]) %>% 
      ungroup() %>% dplyr::select(-state_id, -cat_id, -group_id)
    true_sep = true_data %>% filter(state_id == state_vec[i], cat_id == cat_vec[j]) %>% 
      ungroup() %>% dplyr::select(-state_id, -cat_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    u_sep_methods[,variable_idx,] = 
      getSepUMethods(idx$L6, methods, n = n_sub, f_sep, u_sep, true_sep)
  }
}

display_mat["L6",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))

# L7 ----------------------------------------------------------------------
f_data = pred_L7_all
true_data = df_true_L7

scale2 = unlist(scale2_L7[,3])
weight_vec = (df_weight %>% filter(Level_id == "Level7"))$weight 

f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2

n_state = length(state_vec)
n_dept = length(dept_vec)

true_data_pool = true_data %>% ungroup() %>%dplyr::select(-state_id, -dept_id)
f_pool = f_data %>% arrange(state_id, dept_id, group_id) %>% ungroup() %>% 
  dplyr::select(-state_id, - dept_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

pool_idx = rep(idx$L7, n_state * n_dept) + 
  rep(0:(n_state * n_dept - 1) * n_experts, each = n_sub)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_state * n_dept, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_state * n_dept, n_periods),
                      dimnames = list(methods))

for(i in 1:n_state){
  for(j in 1:n_dept){
    variable_idx = (i-1)*n_dept + j
    f_sep = f_data %>% filter(state_id == state_vec[i], dept_id == dept_vec[j]) %>% 
      ungroup() %>% dplyr::select(-state_id, -dept_id, -group_id)
    true_sep = true_data %>% filter(state_id == state_vec[i], dept_id == dept_vec[j]) %>% 
      ungroup() %>% dplyr::select(-state_id, -dept_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    u_sep_methods[,variable_idx,] = 
      getSepUMethods(idx$L7, methods, n = n_sub, f_sep, u_sep, true_sep)
  }
}

display_mat["L7",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))

# L8 ----------------------------------------------------------------------
f_data = pred_L8_all
true_data = df_true_L8

scale2 = unlist(scale2_L8[,3])
weight_vec = (df_weight %>% filter(Level_id == "Level8"))$weight 

f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2

n_store = length(store_vec)
n_cat = length(cat_vec)

true_data_pool = true_data %>% ungroup() %>%dplyr::select(-store_id, -cat_id)
f_pool = f_data %>% arrange(store_id, cat_id, group_id) %>% ungroup() %>% 
  dplyr::select(-store_id, - cat_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

pool_idx = rep(idx$L8, n_store * n_cat) + 
  rep(0:(n_store * n_cat - 1) * n_experts, each = n_sub)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_store * n_cat, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_store * n_cat, n_periods),
                      dimnames = list(methods))

for(i in 1:n_store){
  for(j in 1:n_cat){
    variable_idx = (i-1)*n_cat + j
    f_sep = f_data %>% filter(store_id == store_vec[i], cat_id == cat_vec[j]) %>% 
      ungroup() %>% dplyr::select(-store_id, -cat_id, -group_id)
    true_sep = true_data %>% filter(store_id == store_vec[i], cat_id == cat_vec[j]) %>% 
      ungroup() %>% dplyr::select(-store_id, -cat_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    u_sep_methods[,variable_idx,] = 
      getSepUMethods(idx$L8, methods, n = n_sub, f_sep, u_sep, true_sep)
  }
}

display_mat["L8",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))


# L9 ----------------------------------------------------------------------
f_data = pred_L9_all
true_data = df_true_L9

scale2 = unlist(scale2_L9[,3])
weight_vec = (df_weight %>% filter(Level_id == "Level9"))$weight 

f_data[,-c(1,2,31)] = f_data[,-c(1,2,31)]/scale2
true_data[,-c(1,2)] = true_data[,-c(1,2)]/scale2

n_store = length(store_vec)
n_dept = length(dept_vec)

true_data_pool = true_data %>% ungroup() %>%dplyr::select(-store_id, -dept_id)
f_pool = f_data %>% arrange(store_id, dept_id, group_id) %>% ungroup() %>% 
  dplyr::select(-store_id, - dept_id, -group_id)
u_pool = as.matrix(f_pool - as.matrix(true_data_pool)  %x% rep(1, n_experts))

pool_idx = rep(idx$L9, n_store * n_dept) + 
  rep(0:(n_store * n_dept - 1) * n_experts, each = n_sub)

u_pool_Linear = getPooledU(pool_idx, "Linear", n_store * n_dept, n_sub, 
                           f_pool, u_pool, true_data_pool)

u_sep_methods = array(NA, dim = c(n_methods, n_store * n_dept, n_periods),
                      dimnames = list(methods))

for(i in 1:n_store){
  for(j in 1:n_dept){
    variable_idx = (i-1)*n_dept + j
    f_sep = f_data %>% filter(store_id == store_vec[i], dept_id == dept_vec[j]) %>% 
      ungroup() %>% dplyr::select(-store_id, -dept_id, -group_id)
    true_sep = true_data %>% filter(store_id == store_vec[i], dept_id == dept_vec[j]) %>% 
      ungroup() %>% dplyr::select(-store_id, -dept_id)
    
    u_sep = t(t(f_sep) - as.numeric(true_sep))
    
    u_sep_methods[,variable_idx,] = 
      getSepUMethods(idx$L9, methods, n = n_sub, f_sep, u_sep, true_sep)
  }
}

display_mat["L9",] = c(getWRMSSE(u_pool_Linear, weight_vec = weight_vec),
                       apply(u_sep_methods, 1, getWRMSSE, weight_vec = weight_vec))

# Print out results -------------------------------------------------------
print(display_mat, digits = 3)
print(apply(display_mat, 1, which.min))
