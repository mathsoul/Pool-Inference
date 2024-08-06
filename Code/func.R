get1SimPoolVar = function(cov_u, n_obs, n_experts, n_products){
  df_train = rmvnorm(n_obs, sigma = cov_u)
  cov_lin = linshrink_cov(df_train)
  cov_inv_lin = solve(cov_lin)
  
  w_df_pool = getPoolWeights(cov_inv_lin, n_products, n_experts)
  
  mean(diag(w_df_pool %*% cov_u %*% t(w_df_pool)))
}

getPoolWeights = function(cov_inv_train, n_products, n_experts){
  B = matrix(NA, n_products, n_products)
  G = matrix(NA, n_products, n_products * n_experts)
  # D = matrix(NA, n_products, n_products)
  for(i in 1:n_products){
    row_idx = 1:n_experts + (i - 1) * n_experts
    for(j in 1:n_products){
      col_idx = 1:n_experts + (j - 1) * n_experts
      B[i,j] = sum(cov_inv_train[row_idx, col_idx])
      # D[i,j] = sum(f_next[row_idx] %*% cov_inv_train[row_idx, col_idx])
      G[i, 1:n_experts + (j-1) * n_experts] = colSums(cov_inv_train[row_idx, col_idx]) 
    }
  }
  
  weights = solve(B) %*% G
  return(weights)
}






getWRMSSE4AllMethods = function(series_idx, data_test, f_data, u_data, n_experts, n_sub, n_items, n_days,
                                min_idx, max_idx){
  true_data = data_test %>% filter(V1 %in% paste0("H", series_idx)) %>% arrange(V1) %>%
    dplyr::select(-V1)
  
  f_comb = f_data %>% filter(id %in% paste0("H", series_idx)) %>%
    arrange(id, group_id) %>% dplyr::select(-id, -group_id)
  
  u_comb = u_data %>% filter(id %in% paste0("H", series_idx)) %>%
    arrange(id, group_id) %>% dplyr::select(-id, -group_id)
  
  min_pool_idx = rep(min_idx, n_items) + rep(0:(n_items - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx, n_items) + rep(0:(n_items - 1) * n_experts, each = n_sub)
  
  u_pool_min_Linear = getPooledUM4(min_pool_idx, "Linear", n_items, n_sub,
                                   f_comb, u_comb, true_data)
  
  u_pool_max_Linear = getPooledUM4(max_pool_idx, "Linear", n_items, n_sub,
                                   f_comb, u_comb, true_data)
  
  u_sep_min_sample = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_max_sample = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_min_Linear = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_max_Linear = matrix(NA, nrow = n_items, ncol = n_days)
  u_EW_min = matrix(NA, nrow = n_items, ncol = n_days)
  u_EW_max = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_min_MS = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_max_MS = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_min_Rob = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_max_Rob = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_min_Cor = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_max_Cor = matrix(NA, nrow = n_items, ncol = n_days)
  
  for(k in 1:n_items){
    f_sep = f_data %>% filter(id %in% paste0("H", series_idx[k])) %>% dplyr::select(-id, - group_id)
    u_sep = u_data %>% filter(id %in% paste0("H", series_idx[k])) %>% dplyr::select(-id, - group_id)
    
    true_sep = data_test %>% filter(V1 %in% paste0("H", series_idx[k])) %>%
      dplyr::select(-V1)
    
    u_sep_min_sample[k,] = unlist(getSepUM4(min_idx, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_min_Linear[k,] = unlist(getSepUM4(min_idx, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_min[k,] = unlist(colMeans(u_sep[min_idx,]))
    u_sep_min_MS[k,] = unlist(getSepUM4(min_idx, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_min_Rob[k,] = unlist(getSepUM4(min_idx, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_min_Cor[k,] = unlist(getSepUM4(min_idx, "Cor", n = n_sub, f_sep, u_sep, true_sep))
    
    u_sep_max_sample[k,] = unlist(getSepUM4(max_idx, "Sample", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_max_Linear[k,] = unlist(getSepUM4(max_idx, "Linear", n = n_sub, f_sep, u_sep, true_sep))
    u_EW_max[k,] = unlist(colMeans(u_sep[max_idx,]))
    u_sep_max_MS[k,] = unlist(getSepUM4(max_idx, "MS", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_max_Rob[k,] = unlist(getSepUM4(max_idx, "Rob", n = n_sub, f_sep, u_sep, true_sep))
    u_sep_max_Cor[k,] = unlist(getSepUM4(max_idx, "Cor", n = n_sub, f_sep, u_sep, true_sep))
  }
  
  display_vec = c(getWRMSSE(u_pool_min_Linear), getWRMSSE(u_sep_min_sample),
                  getWRMSSE(u_sep_min_Linear), getWRMSSE(u_EW_min),
                  getWRMSSE(u_sep_min_MS), getWRMSSE(u_sep_min_Rob),
                  getWRMSSE(u_sep_min_Cor),
                  getWRMSSE(u_pool_max_Linear), getWRMSSE(u_sep_max_sample),
                  getWRMSSE(u_sep_max_Linear), getWRMSSE(u_EW_max),
                  getWRMSSE(u_sep_max_MS), getWRMSSE(u_sep_max_Rob),
                  getWRMSSE(u_sep_max_Cor))
  names(display_vec) = paste(c("Pool Linear","Sep Sample", "Sep Linear", "EW", "Sep MS", "Sep Rob", "Sep Cor"),
                             rep(c("min", "max"), each = 7))
  return(display_vec)
}

getWRMSSE4VarWeighted = function(series_idx, data_test, f_data, u_data, n_experts, n_sub, n_items, n_days,
                                 min_idx, max_idx){
  true_data = data_test %>% filter(V1 %in% paste0("H", series_idx)) %>% arrange(V1) %>%
    dplyr::select(-V1)
  
  f_comb = f_data %>% filter(id %in% paste0("H", series_idx)) %>%
    arrange(id, group_id) %>% dplyr::select(-id, -group_id)
  
  u_comb = u_data %>% filter(id %in% paste0("H", series_idx)) %>%
    arrange(id, group_id) %>% dplyr::select(-id, -group_id)
  
  min_pool_idx = rep(min_idx, n_items) + rep(0:(n_items - 1) * n_experts, each = n_sub)
  max_pool_idx = rep(max_idx, n_items) + rep(0:(n_items - 1) * n_experts, each = n_sub)
  
  u_sep_min_Var = matrix(NA, nrow = n_items, ncol = n_days)
  u_sep_max_Var = matrix(NA, nrow = n_items, ncol = n_days)
  
  for(k in 1:n_items){
    f_sep = f_data %>% filter(id %in% paste0("H", series_idx[k])) %>% dplyr::select(-id, - group_id)
    u_sep = u_data %>% filter(id %in% paste0("H", series_idx[k])) %>% dplyr::select(-id, - group_id)
    
    true_sep = data_test %>% filter(V1 %in% paste0("H", series_idx[k])) %>%
      dplyr::select(-V1)
    
    u_sep_min_Var[k,] = unlist(getSepUM4(min_idx, "Var", n = n_sub, f_sep, u_sep, true_sep))
    
    u_sep_max_Var[k,] = unlist(getSepUM4(max_idx, "Var", n = n_sub, f_sep, u_sep, true_sep))
  }
  
  display_vec = c(getWRMSSE(u_sep_min_Var),getWRMSSE(u_sep_max_Var))
  names(display_vec) = paste(c("Sep Var"),
                             rep(c("min", "max"), each = 1))
  return(display_vec)
}


getMultMinDivIdx = function(idx, items_vec, n_sample, n_sub){
  min_idx = list()
  
  for(i in 1:length(items_vec)){
    min_idx[[i]] = matrix(1:n_sub, n_sample, n_sub, byrow = TRUE)
  }
  return(min_idx)
}

getMultMaxDivIdx = function(u_data, idx, items_vec, n_sample, n_sub){
  max_idx = list()
  
  for(i in 1:length(items_vec)){
    max_idx_1mat = matrix(NA, n_sample, n_sub)
    
    for(j in 1:n_sample){
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

generateMultIdx = function(starting_date, chosen1start, n_samples, items_vec){
  idx = list()
  
  for(j in 1:length(items_vec)){
    idx[[j]] = generate1Idx(starting_date, chosen1start, n_samples, items_vec[j])
  }
  
  return(idx)
}

generate1Idx = function(starting_date, chosen1start, n_samples, n_items){
  idx_candidates = which(starting_date == chosen1start)
  
  map_dfr(1:n_sample,~
            data.frame(matrix(sort(sample(idx_candidates, n_items, replace = FALSE)), nrow = 1)))
}


getCVIdx = function(u_mat, lambda_vec, eta_vec){
  CV_rets = map_dfc(1:length(lambda_vec), ~ getLeaveOneOutRet(u_mat, lambda_vec[.x], eta_vec[.x]))
  
  which.min(apply(CV_rets, 2, function(x) mean(x^2)))
}


getLeaveOneOutRet = function(u_mat, lambda_time_vec, eta_time_vec, sol_mat = NULL){
  n_experts = nrow(u_mat)
  n_days = ncol(u_mat)
  
  if(is.null(sol_mat)){
    sol_mat = getLeaveOneOutSolMat(u_mat, lambda_time_vec, eta_time_vec)
  }
  
  ret_vec = map_dbl(1:n_days, ~
                      getOOSOnePeriodRet(u_mat[,-.x], u_mat[,.x], sol_mat[,.x]))
  
  return(ret_vec)
}

getLeaveOneOutSolMat = function(u_mat, lambda_time_vec, eta_time_vec){
  n_experts = nrow(u_mat)
  n_days = ncol(u_mat)
  
  if(length(lambda_time_vec) == 1){
    lambda_time_vec = rep(lambda_time_vec,n_days)
  }
  
  if(length(eta_time_vec) == 1){
    eta_time_vec = rep(eta_time_vec,n_days)
  }
  
  sol_mat = as.matrix(map_dfc(1:n_days, ~
                                getOneRegSol(u_mat[,-.x], lambda_time_vec[.x], eta_time_vec[.x])))
  
  return(sol_mat)
}

getOOSOnePeriodRet = function(train, test, sol){
  bias_vec = rowMeans(train)
  
  (test - bias_vec) %*% sol
}

getOneRegSol = function(u_train, lambda, eta){
  n_experts = nrow(u_train)
  
  u_demean = scale(t(u_train), center = TRUE, scale = FALSE)
  
  beta_mat = as.matrix(map_dfc(1:n_experts, ~ getRidgepSol(u_demean[,.x], u_demean[,-.x], lambda, eta, -1)))
  var_eps_vec = map_dbl(1:n_experts, ~ sum((u_demean[,.x] - u_demean[,-.x] %*% beta_mat[,.x])^2))
  
  sol = (colSums(beta_mat) - 1)/var_eps_vec
  
  if(abs(sum(sol)) < 1e-5 ){
    sol = rep(1/n_experts, n_experts)
  }else{
    sol = sol/sum(sol)
  }
  return(sol)
}

nZeros = function(vec, threshold){
  length(which(abs(vec) < threshold))
}

getGLASSOpObj = function(S, Theta, lambda, eta = 0){
  log(base::det(Theta)) - sum(diag(S %*% Theta)) - lambda * sum(abs(Theta)) - eta * sum(abs(colSums(Theta)))
}

getTheta22 = function(S22, theta12, Theta11, lambda, eta){
  Q1 = t(theta12) %*% solve(Theta11) %*% theta12
  Q2 = -sum(theta12)
  
  if(Q1 + 1/(S22 + lambda + eta) >= Q2){
    sol = Q1 + 1/(S22 + lambda + eta)
  }else if(Q1 + 1/(S22 + lambda - eta) <= Q2){
    sol = Q1 + 1/(S22 + lambda - eta)
  }else{
    sol = Q2
  }
  
  return(sol)
}

getQuadSol = function(Theta11_inv, s12, lambda, eta = 0, x0 = 0){
  p = ncol(Theta11_inv)
  
  A = rbind(cbind(diag(p),-diag(p), 0 ),
            cbind(-diag(p),-diag(p), 0),
            c(rep(1,p), rep(0, p), -1),
            c(rep(-1,p), rep(0, p), -1))
  
  Q = matrix(0, 2*p + 1, 2*p + 1)
  Q[1:p, 1:p] = 1/2 * Theta11_inv
  
  obj = c(s12, rep(lambda, p), eta)
  
  lb = c(rep(-Inf, p), rep(0, p + 1))
  sense = rep("<=", 2 * p + 2)
  rhs = c(rep(0, 2*p), -x0, x0)
  
  model <- list(A = A, sense = sense, rhs = rhs, lb = lb, obj = obj, Q = Q)
  
  params = list(OutputFlag = 0)
  
  test = gurobi(model, params)
  
  return(test$x[1:p])
}

getRidgepSol = function(y, X, lambda, eta = 0, x0 = 0){
  p = ncol(X)
  n_obs = nrow(X)
  
  A = rbind(c(rep(1,p), -1),
            c(rep(-1,p), -1))
  
  Q = matrix(0, p + 1, p + 1)
  Q[1:p, 1:p] = t(X) %*% X + lambda * diag(p)
  
  obj = c(-2 * t(X) %*% y, eta)
  
  lb = c(rep(-Inf, p),0)
  sense = rep("<=", 2)
  rhs = c(-x0, x0)
  
  model <- list(A = A, sense = sense, rhs = rhs, lb = lb, obj = obj, Q = Q)
  
  params = list(OutputFlag = 0)
  
  test = gurobi(model, params)
  
  return(test$x[1:p])
}

getLASSOpSol = function(y, X, lambda, eta = 0, x0 = 0){
  p = ncol(X)
  n_obs = nrow(X)
  
  A = rbind(cbind(diag(p),-diag(p), 0 ),
            cbind(-diag(p),-diag(p), 0),
            c(rep(1,p), rep(0, p), -1),
            c(rep(-1,p), rep(0, p), -1))
  
  Q = matrix(0, 2*p + 1, 2*p + 1)
  Q[1:p, 1:p] = t(X) %*% X
  
  obj = c(-2 * t(X) %*% y, rep(lambda, p), eta)
  
  lb = c(rep(-Inf, p), rep(0, p + 1))
  sense = rep("<=", 2 * p + 2)
  rhs = c(rep(0, 2*p), -x0, x0)
  
  model <- list(A = A, sense = sense, rhs = rhs, lb = lb, obj = obj, Q = Q)
  
  params = list(OutputFlag = 0)
  
  test = gurobi(model, params)
  
  return(test$x[1:p])
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

MS_method2 <- function(train_data){
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
  
  #I could speed up the following construction but it should be OK.
  #I do not use this optimization that often
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
  
  model$modelsense <- dir
  model$sense      <- c('=', rep('<=', nrow(model$A) - 1))
  model$vtype      <- c(rep('B', n_experts), rep("C", (n_experts-1)*n_experts/2))
  
  params <- list(OutputFlag=1)
  
  result <- gurobi(model, params)
  
  print(round(result$objval/(n_sub * (n_sub -1)/2),4))
  
  return(result$x[1:n_experts])
}

getDistMat = function(var_comb, weight_vec = NULL, dir = "max"){
  n_experts = nrow(var_comb)
  n_items = ncol(var_comb)
  dist_mat = matrix(NA, n_experts, n_experts)
  
  if(is.null(weight_vec)){
    weight_vec = rep(1/n_items, n_items)
  }
  
  for(i in 1:(n_experts-1)){
    for(j in (i+1):(n_experts)){
      dist_mat[j,i] = sqrt((var_comb[i,] - var_comb[j,])^2 %*% weight_vec)
    }
  }
  
  if(dir == "min"){
    identical_idx = which(dist_mat == 0)
    
    if(length(identical_idx) > 0){
      dist_mat[identical_idx] = 10
    }
  }
  
  return(dist_mat)
}

getUCali = function(pred_qtl_all, true_vals, qtl_vec){
  true_vals = as.numeric(true_vals)
  
  n_qtl = length(qtl_vec)
  n_teams = nrow(pred_qtl_all)/n_qtl
  
  cali_mat = matrix(NA, ncol = n_qtl, nrow = n_teams)
  colnames(cali_mat) = qtl_vec
  
  for(i in 1:n_qtl){
    pred_1qtl_all = pred_qtl_all %>% filter(str_detect(id, qtl_vec[i])) %>%
      dplyr::select(-id, -group_id)
    
    err_1qtl_all = t(t(pred_1qtl_all) - true_vals)
    
    cali_mat[,i] = rowMeans(err_1qtl_all > 0)
  }
  return(cali_mat)
}

getUPinball = function(pred_qtl_all, true_vals, qtl_vec){
  true_vals = as.numeric(true_vals)
  
  n_qtl = length(qtl_vec)
  n_teams = nrow(pred_qtl_all)/n_qtl
  
  pinball_mat = matrix(NA, ncol = n_teams, nrow = n_qtl)
  # colnames(pinball_mat) = qtl_vec
  
  for(i in 1:n_qtl){
    pred_1qtl_all = pred_qtl_all %>% filter(str_detect(id, qtl_vec[i])) %>%
      dplyr::select(-id, -group_id)
    
    err_1qtl_all = -t(t(pred_1qtl_all) - true_vals)
    
    pinball_mat[i,] = apply(err_1qtl_all, 1, newloss, tau = as.numeric(qtl_vec[i]))
  }
  
  
  
  # pinball_df = data.frame(pinball_mat) %>% mutate(group_id = 1:50)
  # 
  # pinball_df_long = pivot_longer(pinball_df, cols = "X0.005":"X0.995", values_to = "Loss",
  #                                names_to = "Quantile")
  # 
  # pinball_df_long = pinball_df_long %>% mutate(Quantile = substr(Quantile, 2, 6))
  
  return(pinball_mat)
}


QtlFromWinkler = function(pred_all, true_vals, qtl_vec){
  n_days = ncol(pred_all)
  n_qtl = length(qtl_vec)
  
  err_mat = t(t(pred_all) - as.numeric(true_vals))
  
  qtl_mat = matrix(NA, n_days, n_qtl)
  
  for(i in 1:n_days){
    err_train = err_mat[,-i]
    cov_mat = cov(t(err_train))
    cov_inv = solve(cov_mat)
    w = as.matrix(rowSums(cov_inv)/sum(cov_inv))
    sigma2 = 1/(t(w) %*% cov_inv %*% w)
    
    bias_vec = rowMeans(err_train)
    
    mu = t(w) %*% as.matrix(pred_all[,i] - bias_vec)
    
    qtl_mat[i,] = c(mu) + c(sqrt(sigma2)) * qnorm(qtl_vec)
  }
  
  return(t(qtl_mat))
}

getCov41Qtl = function(pred_qtl_all, qtl){
  pred_1qtl_all = pred_qtl_all %>% filter(str_detect(id, paste0("_", as.character(qtl),"_"))) %>%
    select(-id,-group_id)
  pred_1qtl_centered = scale(pred_1qtl_all, center = TRUE, scale = FALSE)
  cov_1qt = cov(t(pred_1qtl_centered))
  cov_1qt
}


getWRMSSE = function(u_mat, weight_vec = NULL){
  n_items = nrow(u_mat)
  if(is.null(weight_vec)){
    weight_vec = rep(1/n_items, n_items)
  }
  
  sqrt(rowMeans(u_mat^2)) %*% weight_vec
}

getSepUM4 = function(idx, cov_type, n, f_sep, u_sep, true_sep){
  f_comb_agg = map_dfc(1:48, ~ aggForecast(f_sep[idx,.x], u_sep[idx,],
                                           .x, m = 1, n = n, n_days - 1,
                                           cov_type = cov_type))
  u_comb_agg = f_comb_agg - true_sep
  u_comb_agg
}


getSepU = function(idx, cov_type, n, f_sep, u_sep, true_sep){
  f_comb_agg = map_dfc(1:28, ~ aggForecast(f_sep[idx,.x], u_sep[idx,],
                                           .x, m = 1, n = n, n_days - 1,
                                           cov_type = cov_type))
  u_comb_agg = f_comb_agg - true_sep
  u_comb_agg
}

getPooledUM4 = function(idx, cov_type, m, n, f_comb, u_comb, true_data){
  f_comb_agg = map_dfc(1:48, ~ aggForecast(f_comb[idx,.x], u_comb[idx,],
                                           .x, m = m, n = n, n_days - 1,
                                           cov_type = cov_type))
  
  u_comb_agg = f_comb_agg - true_data
  
  u_comb_agg
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

aggForecast = function(f_vec, u_mat, out_idx, m, n, n_days, cov_type = "Sample"){
  u_in = t(u_mat[,-out_idx])
  mean_in = colMeans(u_in)
  #rewrite svd version
  # u_in_demean = scale(u_in, center = TRUE, scale = FALSE)
  # svd_in = svd(u_in_demean, nu = 0, nv = 1)
  
  # n_PC = which.max(map_dbl(1:8,~BaiNG_measure(.x, svd_in$d^2, m*n, n_days)))
  
  # mean plugin
  # plugin_value = mean(svd_in$d[-1]^2)
  # cov_est = (plugin_value * diag(m*n) + (svd_in$d[1]^2 - plugin_value) * svd_in$v %*% t(svd_in$v))/n_days
  
  if(cov_type == "Ville"){
    X = list()
    for(i in 1:m){
      X = c(X, list(unlist(f_vec)[(1:n) + (i-1)*n]))
    }
    
    MEAN_AU_HET = anova_update(X, MEAN = mean, anova = "HET", lme_method = "REML")
    return(MEAN_AU_HET$AU)
    
    # X = list(c(100, 110, 130, 150, 105),
    #          c(90, 85, 76, 94),
    #          c(86, 72, 85, 100, 99))
    
  }
  
  if(cov_type == "NonLin"){
    cov_est = try(nlshrink_cov(u_in))
    if(class(cov_est)[1] == "try-error"){
      print(paste0(out_idx, "NonLin failes"))
      cov_est = diag(ncol(u_in))
    }else{
      if(kappa(cov_est) > 1e10){
        print(paste0(out_idx, "Cond is too large"))
        cov_est = diag(ncol(u_in))
      }
      
    }
    
    # cov_type = "EW"
  }
  
  #PC+EW plugin
  if(cov_type == "EW"){
    one_vec = rep(1, m*n)
    Proj_EW = one_vec - one_vec %*% svd_in$v[,1] %*% t(svd_in$v[,1])
    plugin_value = sum((u_in_demean %*% t(Proj_EW))^2 )/sum(Proj_EW^2)
    cov_est = (plugin_value * diag(m*n) + (svd_in$d[1]^2 - plugin_value) * svd_in$v %*% t(svd_in$v))/n_days
  }
  
  # Factor such that individual variance is matched
  if(cov_type == "Factor"){
    X_PC1 = u_in_demean %*% svd_in$v[,1,drop = FALSE]
    OLS_result = lm(u_in ~ X_PC1)
    
    beta = OLS_result$coefficients[2,, drop = FALSE]
    var_eps = apply(OLS_result$residuals, 2, var)
    
    cov_est = c(var(X_PC1)) * t(beta) %*% beta + diag(var_eps)
  }
  
  if(cov_type == "Var"){
    cov_est = diag(diag(cov(u_in)))
  }
  
  if(cov_type == "Rob"){
    cov_est = cov.rob(u_in)$cov
  }
  
  # NonLin Shrink
  # cov_est = nlshrink_cov(u_in) #it breaks
  
  # Linear Shrink
  if(cov_type == "Linear"){
    cov_est = linshrink_cov(u_in)
  }
  
  if(cov_type == "Sample"){
    cov_est = cov(u_in)
  }
  
  if(cov_type == "Cor"){
    cov_est = getConstCorCov(u_in)
  }
  
  
  if(cov_type == "MS"){
    w_MS = MS_method2(u_in)
    return(w_MS %*% unlist(f_vec - mean_in))
  }
  
  cov_est = cov_est  + 1e-8 * diag(ncol(u_in))
  
  
  if(m > 1){
    f_agg = getMLEFromGurobi(unlist(f_vec - mean_in), solve(cov_est), m, n)
  }else{
    w = rowSums(solve(cov_est))
    w = w/sum(w)
    f_agg = w %*% unlist(f_vec - mean_in)
  }
  
  
  # print(n_PC)
  
  f_agg
}


BaiNG_measure = function(k, eigen_values, n_stocks, train_length){
  log(sum(eigen_values[-(1:k)])) + k *(n_stocks + train_length)/(n_stocks * train_length) *
    log((n_stocks * train_length)/(n_stocks + train_length))
}

getMLEFromGurobi = function(f_vec, Sigma_inv, m, n){
  n_var = m * n
  n_constr = m * (n-1)
  
  A_1 = cbind(1, -diag(n-1))
  A = diag(m) %x% A_1
  
  sense = rep("=", n_constr)
  rhs = rep(0, n_constr)
  
  lb = rep(-Inf, n_var)
  
  obj = -2 * Sigma_inv %*% f_vec
  
  Q = Sigma_inv
  
  model <- list(A = A, sense = sense, rhs = rhs, lb = lb, obj = obj, Q = Q)
  
  params = list(OutputFlag = 1)
  
  test = gurobi(model, params)
  
  test$x[1 + 0:(m-1) * n]
}


simData = function(n_Sim, n_sample, dist = "Normal"){
  if(dist == "Normal"){
    data = matrix(rnorm(n_Sim * n_sample), nrow = n_Sim)
  }else if(dist == "Laplace"){
    data = matrix(rdexp(n_Sim * n_sample), nrow = n_Sim)
  }else if(dist == "t"){
    data = matrix(rt(n_Sim * n_sample, df = 5), nrow = n_Sim)
  }else if(dist == "logistic"){
    data = matrix(rlogis(n_Sim * n_sample), nrow = n_Sim)
  }
  return(data)
}

data2Alpha = function(data){
  row_max = apply(data, 1, max)
  row_mean = rowMeans(data)
  row_SD = apply(data, 1, sd)
  
  # print(cor(row_mean, row_SD))
  
  max_Var = var(row_max)
  mean_Var = var(row_mean)
  SD_Var = var(row_SD)
  
  
  alpha = pnorm(sqrt((max_Var - mean_Var)/SD_Var))
  return(alpha)
}

data2Coefsigma = function(data, alpha = 0.96, sigma){
  qhat = apply(data, 1, quantile, probs = alpha)
  
  # qhat = apply(data, 1, function(x) sort(x)[ncol(data) - 1])
  mean(qhat)/sigma
}



theoLoss = function(mu, alpha, sigma){
  q_true = qnorm(alpha)
  
  q_std = (q_true + mu)/sigma
  
  loss = ((q_true + mu)*pnorm(q_std) + dnorm(q_std) * sigma) *(1 - alpha) +
    (dnorm(q_std) * sigma - (q_true + mu)*(1-pnorm(q_std))) * alpha
  
  return(loss)
}


estDensity = function(err, qtl_vec){
  n_qtl = length(qtl_vec)
  dist_density = density(err, bw = 0.5)
  qtl_vals = quantile(err, qtl_vec)
  
  qtl_idx = map_dbl(1:n_qtl, ~ getDensityIdx(dist_density$x, qtl_vals[.x]))
  # est_density = map_dbl(1:n_qtl, ~ mean(dist_density$y[qtl_idx[.x] + (-10:10)]))
  est_density = dist_density$y[qtl_idx]
  
  # est_cdf = cumsum(dist_density$y)
  # est_cdf = est_cdf/max(est_cdf)
  
  # qtl_idx2 = map_dbl(1:n_qtl, ~max(which(est_cdf < qtl_vec[.x])))
  # est_density2 = dist_density$y[qtl_idx2]
  
  # plot(dist_density$x, dist_density$y)
  # abline(v = qtl_vals, col = 1)
  # abline(v = dist_density$x[qtl_idx2], col = 2)
  
  est_density
}

getDensityIdx = function(density_x, qtl_val){
  which.min(abs(density_x - qtl_val ))
}

estSampleVar = function(pred_vals, true_vals, qtl_vec, n_sample){
  err = as.matrix(pred_vals - true_vals)
  err = err/sd(err)
  est_density = estDensity(err, qtl_vec)
  
  sample_var = qtl_vec * (1-qtl_vec)/est_density^2/n_sample
  names(sample_var) = qtl_vec
  return(sample_var)
}

estSampleplot = function(pred_vals, true_vals, qtl_vec, n_sample){
  err = as.matrix(pred_vals - true_vals)
  err = err/sd(err)
  est_density = density(err, bw = 0.2)#, bw = 0.3
  
  plot(est_density$x, est_density$y)
  abline(v = err)
}

estParaVar = function(pred_vals, true_vals, qtl_vec, n_sample){
  err = as.matrix(pred_vals - true_vals)
  err = as.vector(err/sd(err))
  skewness(err)
  kurtosis(err)
  
  para_var = (1 + 2 * qnorm(qtl_vec) * skewness(err) + qnorm(qtl_vec)^2 * kurtosis(err))/n_sample
  names(para_var) = qtl_vec
  return(para_var)
}



predPara = function(mu, s, n_sample, alpha){
  # (mu + qt(alpha, df = n_sample - 1) * sqrt(1+1/n_sample) * s)
  mu + qnorm(alpha) * s
}

predParaT = function(mu, s, alpha){
  # (mu + qt(alpha, df = n_sample - 1) * sqrt(1+1/n_sample) * s)
  mu + qt(alpha, 5) * s
}


ratioFnormPlot = function(ratio_mat, qtl_vec){
  theme_slides = theme(text=element_text(size=16))
  
  df_ratio = data.frame(t(ratio_mat))
  
  df4plot_long = pivot_longer(df_ratio, cols = 1:9, names_to = "Alpha", values_to = "Ratio")
  df4plot_long$Alpha = factor(rep(qtl_vec, ncol(ratio_mat)), levels = qtl_vec)
  
  p = ggplot(df4plot_long, aes(Alpha, Ratio)) + geom_boxplot() +
    geom_hline(yintercept = 1, linetype = 2) + ylim(0,3.5) + theme_slides
  
  print(p)
}

getSimpleOLS = function(y_vec, X){
  # beta = cov(y_vec, x_vec)/var(x_vec)
  # if(is.nan(beta)){
  #   beta = 0
  # }
  # alpha = mean(y_vec) - beta * mean(x_vec)
  
  coef = lm(y_vec ~ X)$coefficients
}



getOLSSol = function(X, y, Q){
  if(is.null(Q)){
    Q = diag(nrow(X))
  }
  
  solve(t(X) %*% solve(Q) %*% X, t(X) %*% solve(Q) %*% y)
}

EPinBallLoss = function(tau, sample_tau, n_sample){
  q_pop = qnorm(tau)
  q_sample = qnorm(sample_tau)
  
  mu = - q_sample
  sigma2_sample = 1/n_sample * (sample_tau)*(1- sample_tau)/(dnorm(q_sample)^2)
  
  sigma2 = 1 + sigma2_sample
  sigma = sqrt(sigma2)
  
  # mu = 0
  # sigma = 1
  pred_scaled = (0 - mu)/sigma
  
  mu * (tau - pnorm(pred_scaled)) + sigma * dnorm(pred_scaled)
}

# alignRows = function(true_vals, pred_vals, L11_pattern){
#   true_vals = true_vals %>% mutate(id = paste(state_id, item_id, sep = "_"))
#   ids = pred_vals$id[seq(1, nrow(pred_vals), by = 9)]
#   ids = str_extract(ids, L11_pattern)
#
#   df_ids = data.frame(id = ids, no = 1:length(ids))
#
#   join_df = left_join(df_ids, true_vals, by = "id")
#
#   true_vals = join_df %>% dplyr::select(-no, -item_id, -state_id)
# }


alignRowsL12 = function(true_vals, pred_vals){
  true_vals = true_vals %>% mutate(id = paste(item_id, store_id, sep = "_"))
  ids = paste(rep(true_vals$id, each = 9),
              c("0.005", "0.025", "0.165", "0.250", "0.500", "0.750", "0.835", "0.975", "0.995"),
              "evaluation", sep = "_")
  
  df_ids = data.frame(id = ids)
  
  join_df = left_join(df_ids, pred_vals, by = "id")
  
  join_df
}

alignRowsL12Weight = function(true_vals, weight_L12){
  true_vals = true_vals %>% mutate(id = paste(item_id, store_id, sep = "_"))
  weight_L12 = weight_L12 %>% mutate(id = paste(Agg_Level_1, Agg_Level_2, sep = "_"))
  df_ids = data.frame(id = true_vals$id)
  
  join_df = left_join(df_ids, weight_L12, by = "id")
  
  join_df
}

getL10PredLoss = function(true_vals, pred_vals, qtl_vec, n_qtl){
  n_row = nrow(true_vals)
  
  loss_mat = matrix(NA, nrow = n_qtl, ncol = n_row/n_qtl)
  
  for(i in 1:n_qtl){
    idx = seq(i, n_row, by = n_qtl)
    err_mat = as.matrix(true_vals[idx,-1]) - as.matrix(pred_vals[idx,-1])
    loss_mat[i,] = apply(err_mat, 1, newloss, tau = qtl_vec[i])
  }
  
  loss_mat
}

plot4L2IndiQtl = function(pred_L2_all, pred_qtl_L2_all, df_true_L2, state_vec, state_idx,
                          ylim_up = 750){
  true_vals = as.numeric(df_true_L2 %>% filter(str_detect(state_id, state_vec[state_idx])) %>%
                           dplyr::select(-state_id))
  
  loss_P2U_L2 = data.frame()
  
  for(i in 1:50){
    pred_mean = pred_L2_all %>% filter(state_id == state_vec[state_idx], group_id == i) %>%
      dplyr::select(-state_id, -group_id)
    
    pred_P2U_1group = qtlFromMeanLeave1Out(pred_mean, true_vals, qtl_vec, adj_type = "sample")
    loss_P2U_1group = c(getQtlPredLoss(true_vals, pred_P2U_1group, qtl_vec, n_qtl), i, "P2U")
    
    pred_P2UN_1group = qtlFromMeanLeave1Out(pred_mean, true_vals, qtl_vec, adj_type = "normal")
    loss_P2UN_1group = c(getQtlPredLoss(true_vals, pred_P2UN_1group, qtl_vec, n_qtl), i, "P2U_normal")
    
    loss_P2U_L2 = rbind(loss_P2U_L2, loss_P2U_1group, loss_P2UN_1group)
  }
  
  colnames(loss_P2U_L2) = c(qtl_vec, "team_id", "type")
  
  loss_U_L2 = data.frame()
  
  for(i in 1:50){
    pred_U_1group  = pred_qtl_L2_all %>% filter(str_detect(id, state_vec[state_idx]), group_id == i) %>%
      dplyr::select(-id, -group_id)
    loss_U_1group = c(as.numeric(getQtlPredLoss(true_vals, pred_U_1group, qtl_vec, n_qtl)), i, "U")
    
    loss_U_L2 = rbind(loss_U_L2, loss_U_1group)
  }
  
  colnames(loss_U_L2) = c(qtl_vec, "team_id", "type")
  
  df4plot = rbind(loss_P2U_L2, loss_U_L2)
  
  df4plot_long = pivot_longer(df4plot, cols = 1:9, names_to = "Qtl", values_to = "WSPL")
  
  df4plot_long = df4plot_long %>% mutate(Qtl = factor(df4plot_long$Qtl, levels = qtl_vec),
                                         WSPL = as.numeric(WSPL))
  
  ggplot(df4plot_long, aes(x = Qtl, y = WSPL, fill = type)) + geom_boxplot() +
    ylim(c(NA, ylim_up))
  
}

getAggQtlPred = function(pred_all, type = "mean"){
  pred_long = pred_all %>% pivot_longer(cols = F1:F28, names_to = "Days", values_to = "Values") %>%
    mutate(Days = factor(Days, levels = paste0("F", 1:28)))
  
  agg_pred = pred_long %>%  group_by(id, Days) %>% summarise(mean = mean(Values), median = median(Values))
  result = pivot_wider(agg_pred, id_cols = id, names_from = "Days", values_from = type)
}


avgLossAllMean = function(pred_all, true_vals, qtl_vec, n_qtl, pred_type = "factor", qtl_mod_vec = NULL,
                          adj_type = "sample", indi = FALSE, calibration = FALSE){
  if(is.null(qtl_mod_vec)){
    qtl_mod_vec = qtl_vec
  }
  n_experts = nrow(pred_all)
  n_days = ncol(pred_all)
  true_vals = as.numeric(true_vals)
  
  if(pred_type == "factor"){
    svd_result = svd(scale(pred_all, center = TRUE, scale = FALSE))
    
    # n_PC = which.max(map_dbl(1:8,~BaiNG_measure(.x, svd_result$d^2, n_experts, n_days)))
    # stopifnot(n_PC == 1)
    
    mean_matrix = matrix(colMeans(pred_all), byrow = TRUE,
                         nrow = n_experts, ncol = n_days)
    
    pred_all = mean_matrix + svd_result$u[,1] %*% t(svd_result$v[,1]) * svd_result$d[1]
    
    # svd_result = svd(pred_all)
    
    # pred_all = svd_result$u[,1] %*% t(svd_result$v[,1]) * svd_result$d[1]
  }else if(pred_type == "mean"){
    n_experts = 2
    pred_all = matrix(colMeans(pred_all), byrow = TRUE,
                      nrow = n_experts, ncol = n_days)
  }
  if(calibration == FALSE){
    loss_mat = map_dfc(1:n_experts, ~ getQtlPredLoss(true_vals,
                                                     qtlFromMeanLeave1Out(pred_all[.x,], true_vals, qtl_mod_vec, adj_type),
                                                     qtl_vec, n_qtl))
  }else{
    loss_mat = map_dfc(1:n_experts, ~ getCalibration(true_vals,
                                                     qtlFromMeanLeave1Out(pred_all[.x,], true_vals, qtl_mod_vec, adj_type)))
  }
  
  if(indi == FALSE){
    colMeans(loss_mat)
  }else{
    loss_mat
  }
  
}

avgAccLossAllMean = function(pred_all, true_vals, pred_type = "factor"){
  n_experts = nrow(pred_all)
  n_days = ncol(pred_all)
  true_vals = as.numeric(true_vals)
  
  if(pred_type == "factor"){
    svd_result = svd(scale(pred_all, center = TRUE, scale = FALSE))
    
    # n_PC = which.max(map_dbl(1:8,~BaiNG_measure(.x, svd_result$d^2, n_experts, n_days)))
    # stopifnot(n_PC == 1)
    
    mean_matrix = matrix(colMeans(pred_all), byrow = TRUE,
                         nrow = n_experts, ncol = n_days)
    
    pred_all = mean_matrix + svd_result$u[,1] %*% t(svd_result$v[,1]) * svd_result$d[1]
    
    # svd_result = svd(pred_all)
    #
    # pred_all = svd_result$u[,1] %*% t(svd_result$v[,1]) * svd_result$d[1]
  }else if(pred_type == "mean"){
    n_experts = 2
    pred_all = matrix(colMeans(pred_all), byrow = TRUE,
                      nrow = n_experts, ncol = n_days)
  }
  
  loss_mat = map_dfc(1:n_experts, ~ meanSquareRoot(true_vals - pred_all[.x,]))
  
}

getCalibration = function(true_vals, pred_qtl){
  true_vals = as.numeric(true_vals)
  loss_qtl = map_dbl(1:n_qtl, ~ mean(true_vals < pred_qtl[.x,]))
  loss_qtl
}

getQtlPredLoss = function(true_vals, pred_qtl, qtl_vec, n_qtl){
  true_vals = as.numeric(true_vals)
  loss_qtl = map_dbl(1:n_qtl, ~ newloss(true_vals - pred_qtl[.x,], tau = qtl_vec[.x]))
}


absDev = function(x){
  start_idx = min(which(x != 0))
  mean(abs(diff(x[start_idx:length(x)])))
}

squareRoot = function(x){
  start_idx = min(which(x != 0))
  sqrt(mean((diff(x[start_idx:length(x)]))^2))
}

meanSquareRoot = function(x){
  x = as.numeric(x)
  sqrt(mean(x^2))
}

addInfo = function(df, loc_idx, is_Eval = TRUE){
  df_new = cbind(df, loc_idx)
  df_new = df_new %>% mutate(item_id = str_sub(id, 1, V3 - 1),
                             dept_id = str_sub(id, 1, V2 - 1),
                             cat_id = str_sub(id, 1, V1 -1),
                             store_id = str_sub(id, V3 + 1, V5 - 1),
                             state_id = str_sub(id, V3 + 1, V4 - 1),
                             Eval = str_detect(id, "_evaluation"))
  df = cbind(df, df_new %>% dplyr::select(item_id, dept_id, cat_id, store_id, state_id, Eval))
  df = df %>% filter(Eval == is_Eval) %>% dplyr::select(-Eval)
}

lossFromGaba = function(f_mat, y_vec, qtl_vec, calibration = FALSE){
  qtl_pred = qtlFromGabaLeave1Out(f_mat, y_vec, qtl_vec)
  
  err_mat = as.numeric(y_vec) - qtl_pred
  
  n_qtl = length(qtl_vec)
  
  if(calibration == FALSE){
    map_dbl(1:n_qtl, ~ newloss(err_mat[,.x], qtl_vec[.x]))
  }else{
    getCalibration(y_vec, t(qtl_pred))
  }
}

qtlFromGabaLeave1Out = function(f_mat, y_vec, qtl_vec){
  n_days = length(y_vec)
  
  qtl_pred = t(as.matrix(map_dfc(1:n_days, ~ qtlFromGaba(f_mat[,-.x], y_vec[-.x], f_mat[,.x], qtl_vec))))
  
  qtl_pred
}

qtlFromGaba = function(f_train, y_train, f_test, qtl_vec){
  k = length(y_train)
  error_t = t(f_train) - as.numeric(y_train)
  error_t_centered = scale(error_t, center = TRUE, scale = FALSE)
  cor_matrix <- cor(error_t_centered)
  lower_tri <- cor_matrix[lower.tri(cor_matrix)]
  mean_correlation <- mean(lower_tri)
  
  error_t_mean = colMeans(error_t)
  x_bar = mean(f_test-error_t_mean)
  s = sd(f_test-error_t_mean)
  rho = 2*mean_correlation-1
  
  # agg_qtl_pred = x_bar - qt(1-qtl_vec, k) * sqrt((k-1)/k * ((1+rho)/(1-rho) + 1/k)) * s
  agg_qtl_pred = x_bar - qt(1-qtl_vec, k) * 2 * s
}


qtlFromMeanLeave1Out = function(pred_mean, true_vals, qtl_vec, adj_type){
  n_qtl = length(qtl_vec)
  n_days = length(pred_mean)
  # names(pred_mean) = paste("test", 1:n_days)
  
  err = as.matrix(true_vals - pred_mean)
  
  seq_vec = 1:n_days
  names(seq_vec) = paste0("d_", 1:n_days)
  
  adj_mat = map_dfc(seq_vec, ~ getAdjLeave1Out(err, .x, qtl_vec, adj_type))
  
  qtl_mean = t(as.numeric(pred_mean) + t(adj_mat))
  
  return(qtl_mean)
}

getAdjLeave1Out = function(err, test_idx, qtl_vec, adj_type){
  train_idx = -test_idx
  # train_mean = mean(err[train_idx])
  if(adj_type == "sample"){
    adj_train = quantile(err[train_idx], qtl_vec) #- train_mean
  }else if(adj_type == "normal"){
    mu_err = mean(err[train_idx])
    sd_err = sd(err[train_idx])
    z_score = qnorm(qtl_vec)
    # z_score = qt(qtl_vec, length(err) - 1)
    adj_train = mu_err + sd_err * z_score
  }else if(adj_type == "t"){
    mu_err = mean(err[train_idx])
    sd_err = sd(err[train_idx])
    t_score = qt(qtl_vec, 5)
    adj_train = mu_err + sd_err * t_score
  }else if(adj_type == "gamma"){
    fit_gamma = fitdist(unlist(abs(err[train_idx])), distr = "gamma", method = "mme")
    
    qtl_vec_sub = qtl_vec[which(qtl_vec > 0.5)]
    
    adj_sub = qgamma(0.5 + qtl_vec_sub/2, fit_gamma$estimate[1], fit_gamma$estimate[2])
    
    adj_train = c(-rev(adj_sub), 0, adj_sub)
  }else{
    fit_GN = paramp(err[train_idx])
    
    adj_train = qnormp(pr = qtl_vec, mu = fit_GN$mp, sigmap = fit_GN$sp, p = fit_GN$p)
  }
  
  return(adj_train)
}


qtlFromMean = function(pred_mean, true_vals, train_idx, qtl_vec, type = "cheat"){
  n_qtl = length(qtl_vec)
  qtl_mean = matrix(NA, n_qtl, length(pred_mean))
  names(pred_mean) = 1:length(pred_mean)
  
  err = as.matrix(true_vals - pred_mean)
  
  train_result = getMeanSDEst(err[train_idx])
  test_result = getMeanSDEst(err[test_idx])
  if(type == "cheat"){
    adj_train = qnorm(qtl_vec, mean = 0, sd = train_result[2]) - test_result[1]
    adj_test = qnorm(qtl_vec, mean = 0, sd = test_result[2]) - train_result[1]
  }else{
    adj_train = qnorm(qtl_vec, mean = 0, sd = train_result[2]) - train_result[1]
    adj_test = qnorm(qtl_vec, mean = 0, sd = test_result[2]) - test_result[1]
  }
  
  
  qtl_train = map_dfr(1:n_qtl, ~ pred_mean[train_idx] + adj_train[.x])
  qtl_test = map_dfr(1:n_qtl, ~ pred_mean[test_idx] + adj_test[.x])
  
  qtl_mean[,train_idx] = as.matrix(qtl_train)
  qtl_mean[,test_idx] = as.matrix(qtl_test)
  
  return(qtl_mean)
}

newloss = function(err, tau){
  mean(as.numeric((1- tau) * pmax(-err, 0) + tau * pmax(err, 0)))
}

# newlossOld = function(err, tau){
#   (1- tau) * max(-err, 0) + tau * max(err, 0)
# }

# err = as.matrix(df_true_L2[1,-1]) - df_qtl_L2$CA[1,-1]
#
# mean(as.numeric(newloss(err, tau = 0.05)))
#
# mean(newlossOld(err, tau = 0.05))

pullTest = function(pred, mean_pred, train_idx){
  train_dev_mean = mean(as.matrix(pred[train_idx] - mean_pred[train_idx]))
  test_dev_mean = mean(as.matrix(pred[-train_idx] - mean_pred[-train_idx]))
  
  pred[-train_idx] = pred[-train_idx] - (test_dev_mean) + train_dev_mean
  return(pred)
}

maxXinvX = function(x){
  max(x, 1/x)
}

getMeanSDEst = function(x){
  c(median(x), 2*qnorm(0.75)*IQR(x))
  # fitdistr(x, "normal")$estimate
}
