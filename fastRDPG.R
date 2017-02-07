# Fast way of sampling Random dot product graph
# Karl Rohe, Jun Tao, Xintian Han, Norbert Binkiewicz 2017

require(Matrix)
require(igraph)


howManyEdges = function(X,S){
  # this returns the expected number of edges in Poisson gRDPG(X,S) multi-graph.
  # this quantity is proportional to the running time of fastRG
  # Also returns the expected degree and the expected edge density.
  C = diag(colSums(X))
  em = sum(C%*%S%*%C)
  avDeg= em/nrow(X)
  density = em/(nrow(X)^2)
  return(c(em,avDeg, density))
  
}


dcOverlapping = function(theta,pi, B, avgDeg=10){
  # samples from Overlapping SBM with degree correction (DC)
  #  to remove DC, set theta=rep(1,n)
  #  xi  = theta_i * z_i, where [z_i]_j ~ bernoulli(pi_j)
  #  lambda_ij = xi' B xj
  # probability of i connecting to j:  1 - exp(-lambda_ij)
  n = length(theta)
  K = length(pi) 
  
  if(K != nrow(B) | ncol(B) != nrow(B)){
    print("Both dimensions of B must match length of pi")
    return(NA)
  }
  
  
  X = matrix(0, nrow = n, ncol=K)
  for(i in 1:K) X[,i] = rbinom(n,1,pi[i])
  Theta = Diagonal(n,theta)
  X = Theta%*%X
  if(length(avgDeg) == 0) return(fastRG(X,B, simple = T))
  
  em = howManyEdges(X,B)[1]  # expected number of edges in the poisson graph.
  
  scaleIt = em/(n*avgDeg)
  return(fastRG(X,B/scaleIt, simple = T))
}



dcMixed = function(theta,alpha, B, avgDeg=10){
  # samples from mixedmembership SBM with degree correction (DC)
  #  to remove DC, set theta=rep(1,n)
  #  xi  = theta_i * z_i, where z_i ~ dirichlet(alpha)
  #  lambda_ij = xi' B xj
  # probability of i connecting to j:  1 - exp(-lambda_ij)
  n = length(theta)
  K = length(alpha) 
  
  if(K != nrow(B) | ncol(B) != nrow(B)){
    print("Both dimensions of B must match length of alpha")
    return(NA)
  }
  
  
  X = t(sample_dirichlet(n, alpha))
  Theta = Diagonal(n,theta)
  X = Theta%*%X
  
  if(length(avgDeg) == 0) return(fastRG(X,B, simple = T))
  
  em = howManyEdges(X,B)[1]  # expected number of edges in the poisson graph.
  scaleIt = em/(n*avgDeg)
  return(fastRG(X,B/scaleIt, simple = T))
}

dcsbm = function(theta,pi, B, avgDeg=10){
  # samples from a degree corrected stochastic blockmodel
  #  pi a K vector, contains block sampling proportions
  #  theta an n vector, contains degree parameters
  #  B is K time K.
  #   Define lambda_ij = theta[i]*theta[j]*B_{U,V}
  #   where U,V are blockmemberships of i and j, sampled from multinomial(pi)
  
  #   i connects to j with probability 1- exp( -lambda_ij) ,
  #   if lambda_ij is small, then 1- exp( -lambda_ij) \approx lambda_ij
  
  # if avgDeg is set, then B is scaled so that the expected average degree is avgDeg.
  #   in fact, avgDeg is a slight upper bound that is good when the graph is sparse.
  
  #  to make function easy to parameterize, this function does not allow for different 
  #   degree distributions between blocks.
  
  
  n = length(theta)
  K = length(pi) 
  pi = sort(pi/sum(pi))
  
  
  if(K != nrow(B) | ncol(B) != nrow(B)){
    print("Both dimensions of B must match length of pi")
    return(NA)
  }
  
  z = sample(K,n,replace = T, prob = pi)
  # you might want to comment this next line out... but it is here so that pictures are pretty before clustering:
  z = sort(z)  
  X = sparse.model.matrix(~as.factor(z)-1)
  ct= c(0,cumsum(table(z)))
  for(i in 1:K){
    theta[(ct[i]+1):ct[i+1]] = -sort(-theta[(ct[i]+1):ct[i+1]])
  }
  X@x = theta
  
  if(length(avgDeg) == 0) return(fastRG(X,B, simple = T))
  
  em = howManyEdges(X,B)[1]  # expected number of edges in the poisson graph.
  scaleIt = em/(n*avgDeg)
  return(fastRG(X,B/scaleIt, simple = T))
  
}



sbm = function(n,pi, B, avgDeg=NULL){
  # samples from a stochastic blockmodel
  #  pi contains block sampling proportions
  #  B contains probabilities.
  # if avgDeg is set, then B is scaled so that the expected average degree is avgDeg.
  
  
  K = length(pi) 
  # blocks.  
  if(K != nrow(B) | ncol(B) != nrow(B)){
    print("Both dimensions of B must match length of pi")
    return(NA)
  }
  
  z = sample(K,n,replace = T, prob = pi)
  # you might want to comment this next line out... but it is here so that pictures are pretty before clustering:
  z = sort(z)  
  X = sparse.model.matrix(~as.factor(z)-1)
  # there are more arguments that could be specified in fastRG:
  # as set, it will return a simple graph
  #  defaults:  simple = NULL, selfLoops = FALSE, directed = FALSE, multiEdges = TRUE
  S = -log(1-B)  # this is so that bij = 1-exp(-sij)
  
  if(length(avgDeg)==0) return(fastRG(X,S, simple = T))  
  pi = pi/sum(pi)
  
  em = howManyEdges(X,B)[1]  # expected number of edges in the poisson graph.
  
  scaleIt = em/(n*avgDeg)
  S = -log(1-B/scaleIt)
  return(fastRG(X,S, simple = T))
}


fastRG <- function(X, S, 
                   simple = NULL, selfLoops = FALSE, directed = FALSE, multiEdges = TRUE){
  # X is an n x K matrix
  # S is a K x K matrix
  
  # Both X and S should have non-negative entries.
  if(sum(X<0) + sum(S<0)) return(NA)
  
  # if simple = T, then no selfLoops, undirected, and no multiple edges.
  if(length(simple)>0) if(simple){
    selfLoops = FALSE
    directed= FALSE
    multiEdges = FALSE
  }
  
  if(!directed) S = S/2
  # selfLoops = T retains the selfloops, 
  # directed = T does not symmetrize the graph, 
  # multiedges = T keeps the multiple edges.
  
  
  n = nrow(X) 
  K = ncol(X)
  
  C = diag(colSums(X))
  Xt  = apply(X,2,function(x) return(x/sum(x)))
  St = C%*%S%*%C
  m = rpois(1, sum(St)) # number of sampling edges
  if (m==0) return(sparseMatrix(c(1:n),c(1:n),x = 0, dims = c(n, n)))
  
  # to sample U,V from St, we need to sample from a vector, then convert back to matrix.
  UV = sample(K^2, size = m, replace = TRUE, prob = St)  
  tabUV = table(c(1:K^2, UV))
  Uindex = 1:K^2 %% K
  Uindex[Uindex==0]= K
  Vindex = ((1:K^2) %/% K ) + 1
  Vindex = c(1, Vindex)[-(K^2 +1)]
  
  
  #     col_table = table(U)
  #     col_set = as.numeric(as.character(cbind.data.frame(col_table)$U)) 
  #     col_num = cbind.data.frame(col_table)[2]
  #     
  EdgeOut = rep(0, m)
  EdgeIn = EdgeOut
  ticker = 1
  
  for (uv in 1:K^2){
    edgesInThisBlock = tabUV[uv]
    I = sample(n, size = edgesInThisBlock, replace = T, prob = Xt[, Uindex[uv]])
    J = sample(n, size = edgesInThisBlock, replace = T, prob = Xt[, Vindex[uv]])
    
    EdgeOut[ticker:(ticker+edgesInThisBlock-1)] = I
    EdgeIn[ticker:(ticker+edgesInThisBlock-1)] = J
    ticker = ticker + edgesInThisBlock
  }
  
  A = sparseMatrix(EdgeIn, EdgeOut, x = 1, dims = c(n, n))
  if (selfLoops == FALSE){
    diag(A) = 0
  }
  if(!directed) A = A + t(A) 
  if(!multiEdges) A@x[A@x>1]=1
  
  return(A)
}


##### element-wise sampling method
eRDPG <- function(X, dia = TRUE){
  A = matrix(rpois(nrow(X)^2,X%*%t(X)), nrow(X))
  if (dia == FALSE){
    diag(A) = 0
  }
  return(A)
}




