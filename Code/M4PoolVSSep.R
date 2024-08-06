# load("M4 Dataset/ScaledData.Rdata")

# n_sample = 5
# n_sub = 2

# load(paste0("M4 Dataset/SeriesIdx", n_sample,"_",items_vec,".Rdata"))
load(paste0("M4 Dataset/SeriesIdx", n_sample,".Rdata"))

# load(paste0("M4 Dataset/TeamIdx",n_sample, "_", n_sub,"_",items_vec,".Rdata"))
load(paste0("M4 Dataset/TeamIdx",n_sample, "_", n_sub,".Rdata"))
# n_experts = 17
# items_vec = c(3, 7, 9, 10)

# load(paste0("ResultPool/M5_",n_sub,".Rdata"))
# rm(list = c("display_mat_L7","display_mat_L8","display_mat_L9"))

start_time = Sys.time()
# First Starting Date -----------------------------------------------------
WRMSSE_list = foreach(i = 1:length(items_vec),.options.RNG = 20240422,.packages = pkg_names)%dorng%{
  n_items = items_vec[i]
  WRMSSE_1mat = matrix(NA, nrow = n_sample, ncol = 14)
  colnames(WRMSSE_1mat) = paste(c("Pool Linear","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor"),
                                rep(c("min", "max"), each = 7))
  for(j in 1:n_sample){
    series_idx = idx1[[i]][j,]

    WRMSSE_1mat[j, ] = getWRMSSE4AllMethods(series_idx, data_test, f_data, u_data,
                                            n_experts, n_sub, n_items, n_days,
                                            min_idx1[[i]][j,], max_idx1[[i]][j,])
  }
  WRMSSE_1mat
}

save(list = "WRMSSE_list", file = paste0("ResultPool/M4_", n_sample, "_", n_sub,".Rdata"))

# for(k in 1:10){
#   result = foreach(i = 1:10 + (k-1)*10,.options.RNG = 20240422,.packages = pkg_names)%dorng%{
#     series_idx = idx1[[1]][i,]
#     
#     getWRMSSE4AllMethods(series_idx, data_test, f_data, u_data,
#                          n_experts, n_sub, length(series_idx), n_days,
#                          min_idx1[[1]][i,], max_idx1[[1]][i,])
#   }
#   
#   WRMSSE_1mat = matrix(unlist(result), nrow = n_sample/10, byrow = TRUE)
#   colnames(WRMSSE_1mat) = paste(c("Pool Linear","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor"),
#                                 rep(c("min", "max"), each = 7))
#   end_time = Sys.time()
#   print(end_time - start_time)
#   
#   save(list = "WRMSSE_1mat", file = paste0("ResultPool/M4_", n_sample, "_", n_sub, "_",items_vec,"_",k,".Rdata"))
# }

# i= 4
# n_items = items_vec[i]
# for(j in 1:n_sample){
#   print(j)
#   series_idx = idx1[[i]][j,]
# 
#   print(getWRMSSE4Debug(series_idx, data_test, f_data, u_data,
#                   n_experts, n_sub, n_items, n_days,
#                   min_idx1[[i]][j,], max_idx1[[i]][j,]))
# }


# WRMSSE_mat = NULL
# 
# for(k in 1:10){
#   load(paste0("ResultPool/M4_100_3_70_",k,".Rdata"))
#   WRMSSE_mat = rbind(WRMSSE_mat, WRMSSE_1mat)
# }
# 
# WRMSSE_1mat = WRMSSE_mat
# 
# n_sample = 100
# n_sub = 3
# items_vec = 70
# 
# save(list = "WRMSSE_1mat", file = paste0("ResultPool/M4_", n_sample, "_", n_sub, "_",items_vec,".Rdata"))

load("ResultPool/M4_100_2 copy.Rdata")
#
print(c(table(apply(WRMSSE_1mat[,1:7], 1, which.min)),table(apply(WRMSSE_1mat[,8:14], 1, which.min))))
# 
# print(round(colMeans(WRMSSE_1mat), 3))

print(lapply(X = WRMSSE_list, FUN = function(x) c(table(apply(x[,1:7], 1, which.min)),
                                            table(apply(x[,8:14], 1, which.min)))))

xtable(matrix(unlist(lapply(X = WRMSSE_list, FUN = function(x) round(colMeans(x), 3)[c(1,4,2,3,7,5,6) + 7])), ncol = 7,
       byrow = TRUE), digits = 3)

# library(xtable)
# lapply(X = WRMSSE_list, FUN = function(x)
#   xtable(matrix(colMeans(x), byrow = TRUE, nrow =2)[c(2,1),c(1,2,3,7,5,6,4)], digits = 3))
#
# colMeans((WRMSSE_list[[1]]))

xtable(matrix(colMeans(WRMSSE_1mat)[c(1,4,2,3,7,5,6) + 7], nrow = 1), digits = 3)

# cov_est


     