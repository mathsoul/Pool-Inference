# library(tidyverse)
# library(gurobi)

# source("Code/func.R")
# 
# n_sample = 5
# items_vec = c(3, 7, 9, 10)
# n_experts = 17
# n_sub = 3

# load("M4 Dataset/ScaledData.Rdata")
load(paste0("M4 Dataset/SeriesIdx", n_sample, "_", items_vec,".Rdata"))


#First Starting Date ----------------------------------------------------
max_idx1 = getMultMaxDivIdx(u_data, idx1, items_vec, n_sample, n_sub)
min_idx1 = getMultMinDivIdx(idx1, items_vec, n_sample, n_sub)

#Second Starting Date ----------------------------------------------------
# max_idx2 = getMultMaxDivIdx(u_data, idx2, items_vec, n_sample, n_sub)
# min_idx2 = getMultMinDivIdx(idx2, items_vec, n_sample, n_sub)

#Third Starting Date ----------------------------------------------------
# max_idx3 = getMultMaxDivIdx(u_data, idx3, items_vec, n_sample, n_sub)
# min_idx3 = getMultMinDivIdx(idx3, items_vec, n_sample, n_sub)

#Fourth Starting Date ----------------------------------------------------
# max_idx4 = getMultMaxDivIdx(u_data, idx4, items_vec, n_sample, n_sub)
# min_idx4 = getMultMinDivIdx(idx4, items_vec, n_sample, n_sub)

save(list = c(paste0("max_idx",1), paste0("min_idx",1)), 
     file = paste0("M4 Dataset/TeamIdx", n_sample,"_",n_sub,"_",items_vec,".Rdata"))
