sumEachProd = function(u_vec, n_products, n_experts){
  #speed up large matrix multiplication involving matrix E (Proposition 3)
  map_dbl(1:n_products,~sum(u_vec[1:n_experts + (.x-1) * n_experts ]))
}

getLambdaFromLinShrink = function(svd_X, X){ 
  # only use this function when n_products * n_experts >> n_obs
  # it avoids calculating the sample covariance matrix, which could be massive when n_products is large
  # We rescale linshrink_cov to be (n_obs-1) S + lambda I 
  
  n = length(svd_X$d)
  p = nrow(svd_X$v)
  
  m = sum(svd_X$d^2)/(n-1)/p
  
  d2 = (sum(svd_X$d^4)/(n-1)^2 - 2 * m * sum(svd_X$d^2)/(n-1) + m^2 *p)/p
  
  term1 = sum(apply(X, 1, function(x) sum(x^2)^2))
  
  term2 = sum(apply(X, 1, function(x) sum(svd_X$d[-n]^2 * (x %*% svd_X$v)^2)))/(n-1)
  
  term3 = n * sum(svd_X$d^4)/(n-1)^2
  
  b_bar2 = (term1 - 2* term2 + term3)/(n-1)^2/p
  
  b2 = min(d2, b_bar2)
  a2 = d2 - b2
  
  return(b2 * m * (n - 1)/a2)
}


getFastPoolLinearTestU = function(X, X_test, n_products, n_experts, n_obs){
  # We use the Woodbury inequality to derive this algorithm
  # We have checked that it returns the same result as using getPoolWeights
  # X_test has been debiased using train average (mu_vec or mean_in)
  X = scale(X, scale = FALSE) #demean
  svd_X = svd(X, nu = 0, nv = n_obs - 1)
  
  lin_lambda = getLambdaFromLinShrink(svd_X, X)
  
  EtV = apply(svd_X$v, 2, sumEachProd, n_products = n_products, n_experts = n_experts)
  EtU = apply(t(X_test), 2, sumEachProd, n_products = n_products, n_experts = n_experts)
  VtU = t(svd_X$v) %*% t(X_test)
  inv_EtinvE = 1/n_experts *
    (diag(n_products) + EtV %*%
       solve(n_experts * diag(lin_lambda/svd_X$d[-n_obs]^2 + 1) -  t(EtV) %*% EtV) %*% t(EtV))
  
  test = inv_EtinvE %*% (EtU - EtV %*% solve(diag(lin_lambda/svd_X$d[-n_obs]^2 + 1)) %*% VtU)
  
  return(test)
}


simData = function(n_sim, sigma2_x, sigma2_vec, sigma2_eps, n_products, n_experts){
  # When the number of products is large, we could not directly generate data using multivariate normal
  # We have checked that when n_sim is large, the sample covariance matrix is similar to the true covariance matrix
  prod_rand = matrix(rnorm(n_products * n_sim, sd = sqrt(sigma2_x)), nrow = n_sim)
  expert_rand = rmvnorm(n_sim, sigma = diag(sigma2_vec))
  idio_rand = matrix(rnorm(n_sim * n_products * n_experts, sd = sqrt(sigma2_eps)), nrow = n_sim )
  
  data_sim = prod_rand %x% matrix(1, nrow = 1, ncol = n_experts) + 
    matrix(1, nrow = 1, ncol = n_products) %x% expert_rand + idio_rand
  
  return(data_sim)
}

parallelCalM4 = function(j, series_idx_all, team_idx_all, sep_methods,
                         data_test, f_data, u_data, n_periods){
  n_cores = detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Retrieve all function names from the global environment in preparation for parallel computation 
  all_objects <- ls(envir = .GlobalEnv)
  # Filter to include only functions
  all_functions <- all_objects[sapply(mget(all_objects, envir = .GlobalEnv), is.function)]
  
  
  WRMSSE_list = foreach(i = 1:n_sim,.options.RNG = 20240422,
                        .packages = pkg_names, .export = all_functions)%dorng%{
                          series_idx = series_idx_all[[j]][i,]
                          team_idx = team_idx_all[[j]][i,]
                          WRMSSE_vec = getWRMSSE4AllMethodsM4(sep_methods, series_idx, team_idx,
                                                              data_test, f_data, u_data, n_periods)
                          
                          WRMSSE_vec
                        }
  
  WRMSSE_mat = matrix(unlist(WRMSSE_list), nrow = n_sim, byrow = TRUE)
  colnames(WRMSSE_mat) = c("P:Linear", "S:EW", "S:Sample", "S:Linear", "S:Cor",
                           "S:S+EW", "S:Var", "S:Rob")
  
  stopCluster(cl)
  
  WRMSSE_mat
}


getTopIdx = function(idx, prods_vec, n_sim, n_sub){
  min_idx = list()
  
  for(i in 1:length(prods_vec)){
    min_idx[[i]] = matrix(1:n_sub, n_sim, n_sub, byrow = TRUE)
  }
  return(min_idx)
}


OptDivIdx = function(n_experts, n_sub, dist_mat, dir){
  model <- list()
  
  #x1...xn; y11...ynn
  #col:i row:j
  
  # idx_y = n_experts + matrix(1:n_experts^2, n_experts, n_experts)
  idx_y = matrix(NA, n_experts, n_experts)
  count = n_experts
  model$obj <- rep(0, n_experts)
  for(i in 1:(n_experts-1)){
    for(j in (i+1):n_experts){
      count = count + 1
      idx_y[j,i] = count
      
      if(dist_mat[j,i] > 0){
        model$obj = c(model$obj, dist_mat[j,i])
      }else{
        model$obj = c(model$obj, -dist_mat[j,i])
      }
      
    }
  }
  
  
  model$A <- c(rep(1, n_experts), rep(0, (n_experts-1)*n_experts/2))
  model$rhs <- n_sub

  for(i in 1:(n_experts-1)){
    for(j in (i+1):n_experts){
      vec1 = rep(0, n_experts + (n_experts-1)*n_experts/2)
      vec2 = vec1
      vec3 = vec1
      
      vec1[c(i,j)] = 1
      vec1[idx_y[j,i]] = -1
      
      vec2[i] = -1
      vec2[idx_y[j,i]] = 1
      
      vec3[j] = -1
      vec3[idx_y[j,i]] = 1
      
      model$A = rbind(model$A,vec1, vec2, vec3)
      model$rhs = c(model$rhs, c(1,0,0))
    }
  }
  
  model$modelsense <- "max"
  model$sense      <- c('=', rep('<=', nrow(model$A) - 1))
  model$vtype      <- c(rep('B', n_experts), rep("C", (n_experts-1)*n_experts/2))
  
  params <- list(OutputFlag=1)
  
  result <- gurobi(model, params)
  
  print(round(result$objval/(n_sub * (n_sub -1)/2),4))
  
  return(result$x[1:n_experts])
}


getDistMat = function(var_comb, weight_vec = NULL){
  n_experts = nrow(var_comb)
  n_prods = ncol(var_comb)
  dist_mat = matrix(NA, n_experts, n_experts)
  
  if(is.null(weight_vec)){
    weight_vec = rep(1/n_prods, n_prods)
  }
  
  for(i in 1:(n_experts-1)){
    for(j in (i+1):(n_experts)){
      dist_mat[j,i] = sqrt((var_comb[i,] - var_comb[j,])^2 %*% weight_vec)
    }
  }
  
  return(dist_mat)
}


getMultMaxDivIdx = function(u_data, idx, prods_vec, n_sim, n_sub){
  n_experts = 17 #M4 has 17 experts above benchmark
  max_idx = list()
  
  for(i in 1:length(prods_vec)){
    max_idx_1mat = matrix(NA, n_sim, n_sub)
    
    for(j in 1:n_sim){
      series_idx = idx[[i]][j,]
      
      u_comb = u_data %>% filter(id %in% paste0("H", series_idx)) %>%
        arrange(id, group_id) %>% dplyr::select(-id, -group_id)
      
      var_comb = matrix(apply(u_comb, 1, var), nrow = n_experts)
      
      dist_mat_max = getDistMat(var_comb)
      
      max_idx_1mat[j,] = which(round(OptDivIdx(n_experts, n_sub, dist_mat_max, "max")) == 1)
    }
    max_idx[[i]] = max_idx_1mat
  }
  return(max_idx)
}


generate1Idx = function(starting_date, chosen1start, n_sim, n_prods){
  idx_candidates = which(starting_date == chosen1start)
  
  map_dfr(1:n_sim,~
            data.frame(matrix(sort(sample(idx_candidates, n_prods, replace = FALSE)), nrow = 1)))
}


generateMultIdx = function(starting_date, chosen1start, n_sim, prods_vec){
  idx = list()
  
  for(j in 1:length(prods_vec)){
    idx[[j]] = generate1Idx(starting_date, chosen1start, n_sim, prods_vec[j])
  }
  
  return(idx)
}


getWRMSSE4AllMethodsM4 = function(sep_methods, series_idx, team_idx,
                                  data_test, f_data, u_data, n_periods){
  n_prods = length(series_idx)
  n_experts = length(team_idx)
  n_sep_methods = length(sep_methods)
  
  true_data = data_test %>% filter(id %in% paste0("H", series_idx)) %>% arrange(id) %>%
    dplyr::select(-id)
  
  f_comb = f_data %>% filter(id %in% paste0("H", series_idx) & group_id %in% team_idx) %>%
    arrange(id, group_id) %>% dplyr::select(-id, -group_id)
  
  u_comb = u_data %>% filter(id %in% paste0("H", series_idx) & group_id %in% team_idx) %>%
    arrange(id, group_id) %>% dplyr::select(-id, -group_id)
  
  pool_idx = 1:(n_experts * n_prods)
  
  u_pool_Linear = getPooledU(pool_idx, "Linear", n_prods, n_experts, f_comb, u_comb, true_data)
  
  u_sep_methods = array(NA, dim = c(n_sep_methods, n_prods, n_periods), dimnames = list(sep_methods))
  
  for(i in 1:n_prods){
    single_prod_idx = 1:n_experts + (i-1) * n_experts
    
    true_sep = true_data[i,]
    f_sep = tibble(f_comb[single_prod_idx,]) 
    u_sep = u_comb[single_prod_idx,]
    
    u_sep_methods[,i,] = getSepUMethods(1:n_experts, sep_methods, n = n_experts,
                                        f_sep, u_sep, true_sep)
  }
  
  display_vec = c(getWRMSSE(u_pool_Linear), apply(u_sep_methods, 1, getWRMSSE))
  names(display_vec) = c("P:Linear", "S:EW", "S:Sample", "S:Linear", "S:Cor",
                         "S:S+EW", "S:Var", "S:Rob")
  
  display_vec
}

getPoolWeights = function(cov_inv, n_products, n_experts){
  # Formula from Proposition 3 
  E = diag(n_products) %x% matrix(1, n_experts, 1)
  weights = cov_inv %*% E %*% solve(t(E) %*% cov_inv %*% E)
  return(weights)
}

getPooledU = function(idx, cov_type, m, n,f_comb, u_comb, true_data){
  #Leave-one-out result using map_dfc
  f_comb_agg = map_dfc(1:ncol(true_data), ~ aggForecast(f_comb[idx,.x], u_comb[idx,],
                                                        .x, m = m, n = n, n_days - 1,
                                                        cov_type = cov_type))
  
  u_comb_agg = f_comb_agg - true_data
  
  u_comb_agg
}

getSepUMethods = function(idx, methods, n, f_sep, u_sep, true_sep){
  n_methods = length(methods)
  u_comb_mat = matrix(NA, n_methods, ncol(u_sep))
  
  for(i in 1:n_methods){
    f_comb_agg = map_dfc(1:ncol(u_sep), ~ aggForecast(f_sep[idx,.x], u_sep[idx,],
                                                      .x, m = 1, n = n, n_days - 1,
                                                      cov_type = methods[i]))
    
    u_comb_mat[i,] = as.numeric(f_comb_agg - true_sep)
  }
  
  u_comb_mat
} 

get1SimPoolRMSE = function(data_test, n_obs, sigma2_x, sigma2_vec, sigma2_eps, n_products, n_experts){
  # df_train = rmvnorm(n_obs, sigma = cov_u) # computational costly when n_products is large
  data_train = simData(n_obs, sigma2_x, sigma2_vec, sigma2_eps, n_products, n_experts)
  
  mu_vec = colMeans(data_train)
  
  data_test = t(t(data_test) - mu_vec)
  
  if(n_obs > n_products * n_experts){
    cov_est = linshrink_cov(data_train)
    w_pool = getPoolWeights(solve(cov_est), n_products, n_experts)
    u_test = t(data_test %*% w_pool)
  }else{
    u_test = getFastPoolLinearTestU(data_train, data_test, n_products, n_experts, n_obs)
  }

  mean(sqrt(rowMeans(u_test^2)))
}

aggForecast = function(f_vec, u_mat, out_idx, m, n, n_days, cov_type = "Sample"){
  u_in = t(u_mat[,-out_idx])
  mean_in = colMeans(u_in)
  
  if(cov_type == "EW"){
    return(mean(pull(f_vec)))
  }
  
  if(cov_type == "Sample"){
    cov_est = cov(u_in)
  }
  
  if(cov_type == "Linear"){
    cov_est = linshrink_cov(u_in)
  }
  
  if(cov_type == "Cor"){
    cov_est = getConstCorCov(u_in)
  }
  
  if(cov_type == "Var"){
    cov_est = diag(diag(cov(u_in)))
  }
  
  if(cov_type == "S+EW"){
    w_SEW = SoptPlusEW(u_in)
    return(w_SEW %*% unlist(f_vec - mean_in))
  }
  
  if(cov_type == "Rob"){
    cov_est = cov.rob(u_in,seed = 1:20240422)$cov + 1e-8 * diag(n) #to prevent generate a singular covariance estimation
  }
  
  if(m > 1){
    f_agg = t(getPoolWeights(solve(cov_est), m, n)) %*% unlist(f_vec - mean_in)
  }else{
    w = rowSums(solve(cov_est))
    w = w/sum(w)
    f_agg = w %*% unlist(f_vec - mean_in)
  }
  
  f_agg
}

getConstCorCov = function(u_in){
  n = nrow(u_in)
  k = ncol(u_in)
  
  sd_vec = sqrt(diag(cov(u_in)))
  cor_mat = cor(u_in)
  
  lower_tri_idx = which(lower.tri(cor_mat, diag = FALSE))
  
  rho_mean = mean(cor_mat[lower_tri_idx])
  
  rho_est = rho_mean/(1+2*k/(n-1) + 3/(n-1))
  
  const_cor_mat = diag(1-rho_est, k, k)  + rho_est
  
  cov_est = diag(sd_vec) %*% const_cor_mat %*% diag(sd_vec)
  
  return(cov_est)
}

SoptPlusEW <- function(train_data){
  X <- scale(as.matrix(train_data), center = TRUE, scale = FALSE)
  n = nrow(X) - 1
  k = ncol(X)
  ew <- rep(1/ncol(X), ncol(X)) %>% as.matrix()
  l <- matrix(rep(1, ncol(X)), ncol = 1)
  covar <- 1/n*t(X)%*%X + 1e-10 * diag(ncol(train_data))
  covar_inv <- solve(covar)
  ow <- as.numeric((covar_inv%*%l)/(t(l)%*%covar_inv%*%l %>% as.numeric()))
  V <- (n/(n-k+1))*((k-1)/(n-k))*(1/(t(l)%*%covar_inv%*%l))
  BB <- (t(l)%*%covar%*%l)/k^2-n/(n-k+1)*(1/(t(l)%*%covar_inv%*%l))
  B <- max(0, BB)
  lambda <- V/(B+V)
  weight <- c(lambda)*c(ew) + c(1-lambda)*ow
  return(weight)
}

getWRMSSE = function(u_mat, weight_vec = NULL){
  n_prods = nrow(u_mat)
  if(is.null(weight_vec)){
    weight_vec = rep(1/n_prods, n_prods)
  }
  
  sqrt(rowMeans(u_mat^2)) %*% weight_vec
}

extractNSub <- function(input_string) {
  # Use regular expression to extract the number after "Idx"
  match <- regmatches(input_string, regexpr("Idx(\\d+)", input_string))
  number <- as.numeric(sub("Idx", "", match))
  return(number)
}

squareRoot = function(x){
  start_idx = min(which(x != 0))
  sqrt(mean((diff(x[start_idx:length(x)]))^2))
}

# Function to split ID into components
splitId <- function(id) {
  parts <- strsplit(id, "_")[[1]]
  if (length(parts) != 6) {
    stop(paste("ID", id, "does not contain exactly 5 underscores (6 parts)"))
  }
  return(parts)
}
