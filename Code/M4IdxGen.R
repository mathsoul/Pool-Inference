# library(tidyverse)
# source("Code/func.R")

# n_sample = 5


series_info = read.csv("M4 Dataset/M4-info.csv")

n_same_start = sort(table(series_info$StartingDate), decreasing = TRUE)

chosen_idx = which(n_same_start > 10) # L3 has 10 items

chosen_start = names(n_same_start[chosen_idx])


# only 4 names
idx1 = generateMultIdx(series_info$StartingDate, chosen_start[1], n_samples, items_vec)
# idx2 = generateMultIdx(series_info$StartingDate, chosen_start[2], n_samples, items_vec)
# idx3 = generateMultIdx(series_info$StartingDate, chosen_start[3], n_samples, items_vec)
# idx4 = generateMultIdx(series_info$StartingDate, chosen_start[4], n_samples, items_vec)

save(list = c("idx1"), #"idx2" "idx3", "idx4"
     file = paste0("M4 Dataset/SeriesIdx", n_sample,"_", items_vec,".Rdata"))