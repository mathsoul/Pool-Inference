#tobe changed
getPoolWeights = function(cov_inv, n_products, n_experts){
  E = diag(n_products) %x% matrix(1, n_experts, 1)
  weights = cov_inv %*% E %*% solve(t(E) %*% cov_inv %*% E)
  return(weights)
}

getPooledU = function(idx, cov_type, m, n,f_comb, u_comb, true_data){
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

sumEachProd = function(u_vec, n_products, n_experts){
  map_dbl(1:n_products,~sum(u_vec[1:n_experts + (.x-1) * n_experts ]))
}

getLambdaFromLinShrink = function(svd_X, X){ 
  # only use this function when n_products * n_experts > n_obs
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
  prod_rand = matrix(rnorm(n_products * n_sim, sd = sqrt(sigma2_x)), nrow = n_sim)
  expert_rand = rmvnorm(n_sim, sigma = diag(sigma2_vec))
  idio_rand = matrix(rnorm(n_sim * n_products * n_experts, sd = sqrt(sigma2_eps)), nrow = n_sim )
  
  data_sim = prod_rand %x% matrix(1, nrow = 1, ncol = n_experts) + 
    matrix(1, nrow = 1, ncol = n_products) %x% expert_rand + idio_rand
  
  return(data_sim)
}