library(tidyverse)

# Based on the ranking from Submission Info.xlsx
rank_vec = c("118", "245", "237", "072", "069", "036", "078", "260",
             "238", "039", "005", "132", "251", "250", "243", "235", "104")

data_train = read.csv("M4RawData/Hourly-train.csv")
data_test = read.csv("M4RawData/Hourly-test.csv")

stopifnot(all.equal(data_train$V1, paste0("H",1:414)))
stopifnot(all.equal(data_test$V1, paste0("H",1:414)))

# Saled the data according to the scaling used in mean absolute scaled error (MASE)
# documented in [Makridakis et al. (2020)](https://doi.org/10.1016/j.ijforecast.2019.04.014)
n_train = ncol(data_train)
n_period = 24

scale_vec = rowMeans(abs(data_train[,(n_period + 1 + 1):n_train] - 
                           data_train[,(1 + 1):(n_train - n_period)]), na.rm = TRUE)

data_test[,-1] = data_test[,-1]/scale_vec
colnames(data_test) = c("id", paste0("F",1:48))

f_data = data.frame()
u_data = data.frame()

for(i in 1:length(rank_vec)){
  f_mat = read.csv(paste0("M4RawData/submission-", rank_vec[i],".csv"))
  f_mat = f_mat %>% filter(str_detect(id, "H")) 
  
  if(i == 14){
    nums = as.numeric(gsub("\\D", "", f_mat$id)) #The submission file of 250 needs special handling
    f_mat = f_mat[order(nums),]
  }
  
  stopifnot(all.equal(f_mat$id, paste0("H",1:414)))
  
  if(i == 7){
    f_mat = f_mat[,-50] #The submission file of 078 needs special handling
  }
  f_mat[,-1] = f_mat[,-1]/scale_vec
  u_mat = f_mat
  u_mat[,-1] = f_mat[,-1] - data_test[,-1]
  
  f_data = rbind(f_data, data.frame(f_mat, group_id = i))
  u_data = rbind(u_data, data.frame(u_mat, group_id = i))
}


# Select the hourly predictions with same starting time -------------------
series_info = read.csv("M4RawData/M4-info.csv")
series_info = series_info %>% filter(SP == "Hourly") #Only hourly data

n_same_start = sort(table(series_info$StartingDate), decreasing = TRUE)

chosen_start = names(n_same_start[1])

series_info = series_info %>% filter(StartingDate == chosen_start) #With the same start time

u_data = u_data %>% filter(id %in% series_info$M4id) %>% 
  mutate(num = as.numeric(gsub("\\D", "", id))) %>% 
  arrange(num, group_id) %>% dplyr::select(-num)
data_test = data_test %>% filter(id %in% series_info$M4id) 
f_data = f_data %>% filter(id %in% series_info$M4id) %>% 
  mutate(num = as.numeric(gsub("\\D", "", id))) %>%
  arrange(num, group_id) %>% dplyr::select(-num)

save(list = c("u_data", "data_test", "f_data"), file = "CleanedData/M4.Rdata")




