# Fast way of sampling Random dot product graph
# Karl Rohe, Jun Tao, Xintian Han, 2016

rm(list =ls())
library(Matrix)

fastRDPG <- function(X, dia = TRUE){
  
  
  #X, n by K covariate matrix (parameters of nodes in RDPG)
  #dia, Whether to keep self-loops
  #directed,  Whether the edge is directed 
  

  n = nrow(X) #n, The number of nodes in total
  K = ncol(X)

  S = c()
  for (u in 1:K){ S = c(S, sum(X[,u])^2)}
  
  m = rpois(1, sum(S)) # number of sampling edges
  
  U = sample(K, size = m, replace = TRUE, prob = S)
  
  col_table = table(U)
  col_set = as.numeric(as.character(cbind.data.frame(col_table)$U)) 
  col_num = cbind.data.frame(col_table)[2]
  
  Edgein = c()
  Edgeout = c()
  Edgenum = c()
  
  for (col in col_set){
    I = sample(n, size = col_table[[col]], replace = T, prob = X[, col])
    J = sample(n, size = col_table[[col]], replace = T, prob = X[, col])
    Edgein = c(Edgein, I)
    Edgeout = c(Edgeout, J)
    }
  
  pos = n*(Edgein-1) + Edgeout
  pos_table = table(pos)

  pos_set = as.numeric(as.character(cbind.data.frame(pos_table)$pos))

  Edgein = ceiling(pos_set/n)
  Edgeout = pos_set - n*(Edgein-1)
  Edgenum = as.numeric(as.matrix( cbind.data.frame(pos_table)[2]))
    
  A = sparseMatrix(Edgein, Edgeout, x = Edgenum, dims = c(n, n))

  if (dia == FALSE){
    diag(A) = 0
  }
  
  return(A)
}




# ##### naive sampling method
# naiveRDPG <- function(X, dia = TRUE){
#   n = nrow(X)
#   B = X%*%t(X)
#   A = sparseMatrix(i = 1, j = 1, x = 0, dims = c(n,n))
#   for (i in 1:n){
#     for (j in 1:n){
#       A[i,j] = rpois(1, B[i,j])
#     } 
#   }
#   if (dia == FALSE){
#     diag(A) = 0
#   }
#   return(A)
# }
# 
# ########### an easy example
# n = 100
# K = 5
# 
# X = matrix(sample(0:2, size=n*K, replace=T), n, K)
# 
# ptm <- proc.time()
# A = fastRDPG(X, dia = TRUE)
# proc.time() - ptm
# ptm <- proc.time()
# A0 = naiveRDPG(X, dia = TRUE)
# proc.time() - ptm