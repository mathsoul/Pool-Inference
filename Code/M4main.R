source("CodePool/M4LoadPackages.R")
load("M4 Dataset/ScaledData.Rdata")
source("Code/func4Debug.R")

set.seed(20240422)

# items_vec = c(3, 7, 9, 10)
# items_vec = 70
n_experts = 17
n_days = 48
cov_type = "EW"

n_sample = 100
n_sub = 2


for(items_vec in c(3, 7 ,9 ,10)){
  if(items_vec < 70){
    n_core = 10
  }else{
    n_core = 10
  }
  cl <- makeCluster(n_core)
  registerDoParallel(cl)
  
  source("CodePool/M4IdxGen.R")
  # source("CodePool/M4ObtainMaxMinDiversityIdx.R")
  source("CodePool/M4PoolVSSep.R")
  
  stopCluster(cl)
}




