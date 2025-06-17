library(tidyverse)
library(gurobi)

source("Code/func.R")

set.seed(20240422)

prods_vec = c(3, 7, 9, 10, 21, 30, 70)
n_sim = 100
n_experts = 17

series_info = read.csv("M4RawData/M4-info.csv")
series_info = series_info %>% filter(SP == "Hourly") #Only hourly data

n_same_start = sort(table(series_info$StartingDate), decreasing = TRUE)

chosen_start = names(n_same_start[1])

series_info = series_info %>% filter(StartingDate == chosen_start) #With the same start time

series_idx_all = generateMultIdx(series_info$StartingDate, chosen_start, n_sim, prods_vec)

for(n_sub in c(2,3,10,15)){
  top_idx_all = getTopIdx(series_idx_all, prods_vec, n_sim, n_sub)
  max_idx_all = getMultMaxDivIdx(u_data, series_idx_all, prods_vec, n_sim, n_sub)
  
  team_idx_all = top_idx_all
  save(list = c("series_idx_all", "team_idx_all"),
       file = paste0("TeamIndex/M4TopIdx", n_sub, ".Rdata"))
  
  team_idx_all = max_idx_all
  save(list = c("series_idx_all", "team_idx_all"),
       file = paste0("TeamIndex/M4MaxIdx", n_sub, ".Rdata"))
}
