library(tidyverse)
library(purrr)
library(tidyr)
library(ggpubr)

# This file generates Figure 1

cvec = c(1,1,1) 
sigma2_1 = 1
sigma2_2 = 1
sigma2_eps = 10
one_3 = rep(1,3)

# We use Eq. (9) to obtain the optimal weights of pooled inference
getWMatPool = function(sigma2_1, sigma2_2, sigma2_eps, cvec){
  sigma2_exp = c(sigma2_1, sigma2_2)
  Sigma_u = diag(3) %x% matrix(1,2,2) + (cvec %*% t(cvec)) %x% diag(sigma2_exp) +
    sigma2_eps * diag(6)
  E = diag(3) %x% rep(1,2)
  wMat = solve(Sigma_u) %*% E %*% solve(t(E) %*% solve(Sigma_u) %*% E) 
  return(wMat)
}

# We use Eq. (EC.10) to obtain the optimal weights of separate inference
getWMatSep = function(sigma2_1, sigma2_2, sigma2_eps, cvec){
  sigma2_exp = c(sigma2_1, sigma2_2)
  Sigma_u = diag(3) %x% matrix(1,2,2) + (cvec %*% t(cvec)) %x% diag(sigma2_exp) +
    sigma2_eps * diag(6)
  wMat = matrix(0, 6, 3)
  for(i in 1:3){
    idx = 1:2 + (i-1)*2
    wsol = rowSums(solve(Sigma_u[idx, idx]))
    wMat[idx,i] = wsol/sum(wsol)
  }
  return(wMat)
}


sigma2_2_vec = seq(1, 100, by = 0.5)

# 1,3,5 columns are the weights of the first experts for variables 1,2,3
w_mat_prod1_pool = map_dfc(sigma2_2_vec, ~getWMatPool(sigma2_1, .x, sigma2_eps, cvec)[c(1,3,5),1])
w_mat_prod1_sep = map_dfc(sigma2_2_vec, ~getWMatSep(sigma2_1, .x, sigma2_eps, cvec)[c(1,3,5),1])

w_mat_prod2_pool = map_dfc(sigma2_2_vec, ~getWMatPool(sigma2_1, .x, sigma2_eps, cvec)[c(1,3,5),2])
w_mat_prod2_sep = map_dfc(sigma2_2_vec, ~getWMatSep(sigma2_1, .x, sigma2_eps, cvec)[c(1,3,5),2])

w_mat_prod3_pool = map_dfc(sigma2_2_vec, ~getWMatPool(sigma2_1, .x, sigma2_eps, cvec)[c(1,3,5),3])
w_mat_prod3_sep = map_dfc(sigma2_2_vec, ~getWMatSep(sigma2_1, .x, sigma2_eps, cvec)[c(1,3,5),3])


df4plot = data.frame(Sigma2_2 = sigma2_2_vec,
                     Separate = as.numeric(w_mat_prod1_sep[1,]),
                     Pooled = as.numeric(w_mat_prod1_pool[1,]))

df4plot_long = pivot_longer(df4plot, cols = Separate:Pooled, names_to = "Type",
                            values_to = "Weight")

p1 = ggplot(df4plot_long, aes(Sigma2_2, Weight, linetype = Type)) + geom_line(linewidth = 1.25) + 
  ylim(0,1) + xlab(expression(sigma[2]^2)) + ylab(expression(w['1,1'])) +
  theme_bw(18) + theme(legend.text = element_text(size = 25)) + 
  scale_linetype_discrete(name = "")


df4plot = data.frame(Sigma2_2 = sigma2_2_vec,
                     Separate = as.numeric(w_mat_prod1_sep[2,]),
                     Pooled = as.numeric(w_mat_prod1_pool[2,]))

df4plot_long = pivot_longer(df4plot, cols = Separate:Pooled, names_to = "Type",
                            values_to = "Weight")

p2 = ggplot(df4plot_long, aes(Sigma2_2, Weight, linetype = Type)) + geom_line(linewidth = 1.25) + 
  ylim(0,1) + xlab(expression(sigma[2]^2)) + ylab(expression(w['2,1'] == w['3,1'])) +
  theme_bw(18) + theme(legend.text = element_text(size = 25)) + 
  scale_linetype_discrete(name = "")


df4plot = data.frame(Sigma2_2 = sigma2_2_vec,
                     Separate = as.numeric(colSums(w_mat_prod1_sep)),
                     Pooled = as.numeric(colSums(w_mat_prod1_pool)))


df4plot_long = pivot_longer(df4plot, cols = Separate:Pooled, names_to = "Type",
                            values_to = "Weight")


p3 = ggplot(df4plot_long, aes(Sigma2_2, Weight, linetype = Type)) + geom_line(linewidth = 1.25) + 
  ylim(0,1) + xlab(expression(sigma[2]^2)) + ylab(expression(w['1,1'] + w['2,1'] + w['3,1'])) +
  theme_bw(18) + theme(legend.text = element_text(size = 25)) + 
  scale_linetype_discrete(name = "")


# pdf("FigurePool/wMatSigma2Prod1.pdf", width = 16, height = 6)
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="top")
# dev.off()