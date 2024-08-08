#tobe changed
getPoolWeights = function(cov_inv, n_products, n_experts){
  E = diag(n_products) %x% matrix(1, n_experts, 1)
  weights = cov_inv %*% E %*% solve(t(E) %*% cov_inv %*% E)
  return(weights)
}

getPooledU = function(idx, cov_type, m, n,f_comb, u_comb, true_data){
  f_comb_agg = map_dfc(1:28, ~ aggForecast(f_comb[idx,.x], u_comb[idx,],
                                           .x, m = m, n = n, n_days - 1,
                                           cov_type = cov_type))
  if(ncol(true_data) == 29){
    u_comb_agg = f_comb_agg - true_data[,-c(1)]
  }else{
    u_comb_agg = f_comb_agg - true_data[,-c(1,2)]
  }
  
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

get1SimPoolVar = function(cov_u, n_obs, n_products, n_experts){
  df_train = rmvnorm(n_obs, sigma = cov_u)
  cov_lin = linshrink_cov(df_train)
  cov_inv_lin = solve(cov_lin)
  
  w_df_pool = getPoolWeights(cov_inv_lin, n_products, n_experts)
  
  mean(diag(w_df_pool %*% cov_u %*% t(w_df_pool)))
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
  X <- as.matrix(train_data)
  n = nrow(train_data)
  k = ncol(train_data)
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
  n_items = nrow(u_mat)
  if(is.null(weight_vec)){
    weight_vec = rep(1/n_items, n_items)
  }
  
  sqrt(rowMeans(u_mat^2)) %*% weight_vec
}

extractNSub <- function(input_string) {
  # Use regular expression to extract the number after "Idx"
  match <- regmatches(input_string, regexpr("Idx(\\d+)", input_string))
  number <- as.numeric(sub("Idx", "", match))
  return(number)
}
