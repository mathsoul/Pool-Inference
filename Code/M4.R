library(tidyverse)
library(nlshrink)
library(MASS)
library(xtable)
# Parallel computation
library(MCMCpack)
library(doParallel)
library(doRNG)
pkg_names = c("tidyverse", "nlshrink", "MASS")

source("Code/func.R")

load("CleanedData/M4.Rdata")

team_idx_name = "TeamIndex/M4TopIdx2.Rdata" #Table 7. Takes 3 hours using 7 cores.
# team_idx_name = "TeamIndex/M4TopIdx15.Rdata" #Table 8. Takes 5 hours using 7 cores.
# team_idx_name = "TeamIndex/M4TopIdx3.Rdata" #Table EC.3. Takes 4 hours using 7 cores.
# team_idx_name = "TeamIndex/M4TopIdx10.Rdata" #Table EC.4. Takes 4 hours using 7 cores.
# team_idx_name = "TeamIndex/M4MaxIdx2.Rdata" #Table EC.7. Takes 3 hours using 7 cores.
# team_idx_name = "TeamIndex/M4MaxIdx15.Rdata" #Table EC.8. Takes 5 hours using 7 cores.

load(team_idx_name)

n_sim = 100 #number of simulations
n_periods = 48 #number of observations

prods_vec = c(3, 7, 9, 10, 21, 30, 70)

sep_methods = c("EW", "Sample", "Linear", "Cor", "S+EW", "Var", "Rob")

display_mat = matrix(NA, nrow = length(prods_vec), ncol = 1 + length(sep_methods))
win_pct_vec = rep(NA, length(prods_vec))
names(win_pct_vec) = prods_vec

for(j in 1:length(prods_vec)){
  start_time = Sys.time()
  
  # parallel computation might fail for Windows
  WRMSSE_mat = parallelCalM4(j, series_idx_all, team_idx_all, sep_methods,
                             data_test, f_data, u_data, n_periods)
  
  display_mat[j,] = round(colMeans(WRMSSE_mat),3)
  win_pct_vec[j] = table(apply(WRMSSE_mat, 1, which.min))[1]/n_sim
  
  end_time = Sys.time()
  print(end_time - start_time)
}

rownames(display_mat) = prods_vec

print(display_mat, digits = 3)
win_pct_vec




