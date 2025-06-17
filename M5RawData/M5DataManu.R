library(tidyverse)
library(data.table)
source("Code/func.R")

# True Values --------------------------------------------------
df_true = fread("M5RawData/sales_test_evaluation.csv")

df_true_L1 = df_true %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L2 = df_true %>% group_by(state_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L3 = df_true %>% group_by(store_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L4 = df_true %>% group_by(cat_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L5 = df_true %>% group_by(dept_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L6 = df_true %>% group_by(state_id, cat_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L7 = df_true %>% group_by(state_id, dept_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L8 = df_true %>% group_by(store_id, cat_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 
df_true_L9 = df_true %>% group_by(store_id, dept_id) %>% summarise(across(d_1942:d_1969, ~ sum(.x))) 

dept_vec = sort(unique(df_true$dept_id))
cat_vec = sort(unique(df_true$cat_id))
store_vec = sort(unique(df_true$store_id))
state_vec = sort(unique(df_true$state_id))

# Weight Info -------------------------------------------------------------
df_weight = fread("M5RawData/weights_evaluation.csv")
weight_tags = paste0("Level", 1:12, "$", collapse = "|")
df_weight = df_weight %>% filter(str_detect(Level_id, weight_tags))

# Scaling Info ------------------------------------------------------------
df_true_train_all = fread("M5RawData/sales_train_evaluation.csv")

df_train_L1 = df_true_train_all%>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L2 = df_true_train_all%>% group_by(state_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L3 = df_true_train_all%>% group_by(store_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L4 = df_true_train_all%>% group_by(cat_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L5 = df_true_train_all%>% group_by(dept_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L6 = df_true_train_all%>% group_by(state_id, cat_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L7 = df_true_train_all%>% group_by(state_id, dept_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L8 = df_true_train_all%>% group_by(store_id, cat_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))
df_train_L9 = df_true_train_all%>% group_by(store_id, dept_id) %>% summarise(across(d_1:d_1941, ~ sum(.x)))

scale2_L1 = apply(df_train_L1, 1, squareRoot)
scale2_L2 = cbind(df_train_L2[,1], value = apply(df_train_L2[,-1], 1, squareRoot))
scale2_L3 = cbind(df_train_L3[,1], value = apply(df_train_L3[,-1], 1, squareRoot))
scale2_L4 = cbind(df_train_L4[,1], value = apply(df_train_L4[,-1], 1, squareRoot))
scale2_L5 = cbind(df_train_L5[,1], value = apply(df_train_L5[,-1], 1, squareRoot))

scale2_L6 = cbind(df_train_L6[,c(1,2)], value = apply(df_train_L6[,-c(1,2)], 1, squareRoot))
scale2_L7 = cbind(df_train_L7[,c(1,2)], value = apply(df_train_L7[,-c(1,2)], 1, squareRoot))
scale2_L8 = cbind(df_train_L8[,c(1,2)], value = apply(df_train_L8[,-c(1,2)], 1, squareRoot))
scale2_L9 = cbind(df_train_L9[,c(1,2)], value = apply(df_train_L9[,-c(1,2)], 1, squareRoot))


# Read all submissions ----------------------------------------------------
folder_name = "Submissions"

file_names_vec = list.files(paste0("M5RawData//", folder_name, "/"))

n_files = length(file_names_vec)

pred_L1_all = data.frame()
pred_L2_all = data.frame()
pred_L3_all = data.frame()
pred_L4_all = data.frame()
pred_L5_all = data.frame()
pred_L6_all = data.frame()
pred_L7_all = data.frame()
pred_L8_all = data.frame()
pred_L9_all = data.frame()

for(i in 1:n_files){
  
  file_name = paste0("M5RawData/", folder_name, "/", file_names_vec[i])
  df = fread(file_name)
  
  # Apply the function to all IDs and create a data frame
  id_info <- data.frame(do.call(rbind, lapply(df$id, splitId)))
  colnames(id_info) <- c("cat_id", "dept_id", "item_id", "state_id", "store_id", "usage")
  
  id_info = id_info %>% mutate(dept_id = paste0(cat_id,"_", dept_id),
                             store_id = paste0(state_id,"_", store_id))
  
  df = cbind(df, id_info) %>% filter(usage == "evaluation")
  
  df_L1 = df %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L2 = df %>% group_by(state_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L3 = df %>% group_by(store_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L4 = df %>% group_by(cat_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L5 = df %>% group_by(dept_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L6 = df %>% group_by(state_id, cat_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L7 = df %>% group_by(state_id, dept_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L8 = df %>% group_by(store_id, cat_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  df_L9 = df %>% group_by(store_id, dept_id) %>% summarise(across(F1:F28, ~ sum(.x))) %>% mutate(group_id = i)
  
  stopifnot(nrow(df_L2) == 3)
  
  pred_L1_all = rbind(pred_L1_all, df_L1)
  pred_L2_all = rbind(pred_L2_all, df_L2)
  pred_L3_all = rbind(pred_L3_all, df_L3)
  pred_L4_all = rbind(pred_L4_all, df_L4)
  pred_L5_all = rbind(pred_L5_all, df_L5)
  pred_L6_all = rbind(pred_L6_all, df_L6)
  pred_L7_all = rbind(pred_L7_all, df_L7)
  pred_L8_all = rbind(pred_L8_all, df_L8)
  pred_L9_all = rbind(pred_L9_all, df_L9)
}

save(list = c(paste0("df_true_L",1:9), paste0("scale2_L", 1:9), paste0("pred_L", 1:9,"_all"),
              "df_weight","dept_vec", "cat_vec", "store_vec", "state_vec"), 
     file = "CleanedData/M5.Rdata")
