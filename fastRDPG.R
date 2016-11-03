# Fast way of sampling Random dot product graph
# Karl Rohe, Jun Tao, Xintian Han, 2016

rm(list =ls())
library(Matrix)

fRDPG <- function(X, dia = TRUE){
  
  #X, n by K covariate matrix (parameters of nodes in RDPG)
  #dia, Whether to keep self-loops
  #directed,  Whether the edge is directed 
  n = nrow(X) #n, The number of nodes in total
  K = ncol(X)
  
  S = colSums(X)^2
  m = rpois(1, sum(S)) # number of sampling edges
  if (m==0){
    return(sparseMatrix(c(1:n),c(1:n),x = 0, dims = c(n, n)))}
  else{
    U = sample(K, size = m, replace = TRUE, prob = S)
    col_table = table(U)
    col_set = as.numeric(as.character(cbind.data.frame(col_table)$U)) 
    col_num = cbind.data.frame(col_table)[2]
    
    EdgeOut = rep(0, m)
    EdgeIn = EdgeOut
    ticker = 1
    
    for (col in col_set){
      I = sample(n, size = col_table[[col]], replace = T, prob = X[, col])
      J = sample(n, size = col_table[[col]], replace = T, prob = X[, col])
      
      EdgeOut[ticker:(ticker+col_table[[col]]-1)] = I
      EdgeIn[ticker:(ticker+col_table[[col]]-1)] = J
      ticker = ticker + col_table[[col]]
    }
    
    A = sparseMatrix(EdgeIn, EdgeOut, x = 1, dims = c(n, n))
    if (dia == FALSE){
      diag(A) = 0
    }
    
    return(A)
  }
}

##### element-wise sampling method
eRDPG <- function(X, dia = TRUE){
  A = matrix(rpois(nrow(X)^2,X%*%t(X)), nrow(X))
  if (dia == FALSE){
    diag(A) = 0
  }
  return(A)
}