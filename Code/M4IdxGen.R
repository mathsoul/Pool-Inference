library(tidyverse)
library(gurobi)

source("Code/func.R")

load("CleanedData/M4.Rdata")

set.seed(20240422)

prods_vec = c(3, 7, 9, 10, 21, 30, 70)
n_sim = 100
n_sub = 15
n_experts = 17

series_info = read.csv("RawData/M4-info.csv")

n_same_start = sort(table(series_info$StartingDate), decreasing = TRUE)

chosen_start = names(n_same_start[1])

series_idx_all = generateMultIdx(series_info$StartingDate, chosen_start[1], n_sim, prods_vec)

idx = series_idx_all

top_idx_all = getTopIdx(series_idx_all, prods_vec, n_sim, n_sub)
max_idx_all = getMultMaxDivIdx(u_data, series_idx_all, prods_vec, n_sim, n_sub)

team_idx_all = top_idx_all
save(list = c("series_idx_all", "team_idx_all"),
     file = paste0("TeamIndex/M4TopIdx", n_sub, ".Rdata"))

team_idx_all = max_idx_all
save(list = c("series_idx_all", "team_idx_all"),
     file = paste0("TeamIndex/M4MaxIdx", n_sub, ".Rdata"))
