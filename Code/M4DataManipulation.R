library(tidyverse)

rank_vec = c("118", "245", "237", "072", "069", "036", "078", "260",
             "238", "039", "005", "132", "251", "250", "243", "235", "104")
data_train = read.csv("Data/Hourly-train.csv")
data_test = read.csv("Data/Hourly-test.csv")

stopifnot(all.equal(data_train$V1, paste0("H",1:414)))
stopifnot(all.equal(data_test$V1, paste0("H",1:414)))

n_train = ncol(data_train)
n_period = 24

scale_vec = rowMeans(abs(data_train[,(n_period + 1 + 1):n_train] - 
                           data_train[,(1 + 1):(n_train - n_period)]), na.rm = TRUE)

data_test[,-1] = data_test[,-1]/scale_vec

f_data = data.frame()
u_data = data.frame()

for(i in 1:length(rank_vec)){
  f_mat = read.csv(paste0("M4 Dataset/submission-", rank_vec[i],".csv"))
  stopifnot(all.equal(f_mat$id, paste0("H",1:414)))
  f_mat[,-1] = f_mat[,-1]/scale_vec
  u_mat = f_mat
  u_mat[,-1] = f_mat[,-1] - data_test[,-1]
  
  f_data = rbind(f_data, data.frame(f_mat, group_id = i))
  u_data = rbind(u_data, data.frame(u_mat, group_id = i))
}

save(list = c("u_data", "data_test", "f_data"), file = "M4 Dataset/ScaledData.Rdata")




