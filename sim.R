library(ks)
library(abind)
library(meanShiftR)
source("kms_mod.R")

# Data : list of matrix data
# K : number of neighbors
# H : bandwidth matrix
# method : FB (fixed bandwidth), SP (sample point), NN (Nearest neighbors) 

sim_matrix_clust <- function(data, H, K, method = c("FB", "SP", "NN")) {
  
  t <- dim(data)[1]
  p <- dim(data)[2]
  n <- dim(data)[3] 
  x_mat <- t(apply(data, 3, c))
  tol_NNMS <- 0.01 * max(apply(x_mat, 2, IQR))*ncol(x_mat)
  
  if (method == "FB") {
    if (missing(H)) {
      H <- hns(x_mat, deriv.order = 1)*diag(p*t)
    } 
    res <- kms(x_mat, H = H, merge = TRUE, min.clust.size = 10)
  }
  
  if (method == "SP") {
    if (missing(H)) {
      H <- hns(x_mat, deriv.order = 1)*diag(p*t)
    } 
    if (missing(K)) {
      K <- 5*sqrt(n)
    }
    res <- kms.mod(x_mat, H = H, K = K, merge = TRUE, min.clust.size = 10, kernel = "Adaptive")
  }
  
  if (method == "NN") {
    if (missing(K)) {
      K <- 5*sqrt(n)
    }
    res <- NNMS(x_mat, K = K, tol = tol_NNMS)
  }
  
  return(res)
}

