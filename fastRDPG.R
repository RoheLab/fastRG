# Fast way of sampling Random dot product graph
# Karl Rohe, Jun Tao, Xintian Han, Norbert Binkiewicz 2017

require(Matrix)
require(igraph)  # this is needed for the function sample_dirichlet


howManyEdges = function(X,S, Y=NULL){
  # this returns the expected number of edges in Poisson gRDPG(X,S) multi-graph.
  # this quantity is proportional to the running time of fastRG
  # Also returns the expected degree and the expected edge density.
  if(is.null(Y)) Y<-X
  Cx = diag(colSums(X))
  Cy = diag(colSums(Y))
  em = sum(Cx%*%S%*%Cy)
  avDeg= em/nrow(X)
  density = em/(nrow(X)*nrow(Y))
  return(c(em,avDeg, density))
  
}


dcOverlapping = function(theta,pi, B, ...){
  # samples from Overlapping SBM with degree correction (DC)
  #  to remove DC, set theta=rep(1,n)
  #  xi  = theta_i * z_i, where [z_i]_j ~ bernoulli(pi_j)
  #  lambda_ij = xi' B xj
  # probability of i connecting to j:  1 - exp(-lambda_ij)
  if(length(theta)==1){
    n = theta
    theta = rep(1, n)
  }
  if(length(theta)>1) n = length(theta)
  K = length(pi) 
  
  if(K != nrow(B) | ncol(B) != nrow(B)){
    print("Both dimensions of B must match length of pi")
    return(NA)
  }
  
  
  X = matrix(0, nrow = n, ncol=K)
  for(i in 1:K) X[,i] = rbinom(n,1,pi[i])
  X = X[order(X%*%(1:K)),]
  Theta = Diagonal(n,theta)
  X = Theta%*%X

  return(fastRG(X,B, ...))
}



dcMixed = function(theta,alpha, B, ...){
  # samples from mixedmembership SBM with degree correction (DC)
  #  to remove DC, set theta=rep(1,n)
  #  xi  = theta_i * z_i, where z_i ~ dirichlet(alpha)
  #  lambda_ij = xi' B xj
  # probability of i connecting to j:  1 - exp(-lambda_ij)
  if(length(theta)==1){
    n = theta
    theta = rep(1, n)
  }
  if(length(theta)>1) n = length(theta)
  K = length(alpha) 
  
  if(K != nrow(B) | ncol(B) != nrow(B)){
    print("Both dimensions of B must match length of alpha")
    return(NA)
  }
  
  
  X = t(sample_dirichlet(n, alpha))
  X = X[order(X%*%(1:K)),]
  Theta = Diagonal(n,theta)
  X = Theta%*%X
  
  return(fastRG(X,B, ...))
}

dcsbm = function(theta,pi, B, ...){
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
  B = B[order(pi), ]; B = B[, order(pi)]
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
  
  return(fastRG(X,B, ...))
  
}







dcsbm = function(theta,pi, B, ...){
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
  B = B[order(pi), ]; B = B[, order(pi)]
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
  
  return(fastRG(X,B, ...))
  
}










sbm = function(n,pi, B, PoissonEdges = F, avgDeg =NULL, ...){
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
  X = model.matrix(~factor(as.character(z), levels = as.character(1:K))-1)
  # there are more arguments that could be specified in fastRG:
  # as set, it will return a simple graph
  #  defaults:  simple = NULL, selfLoops = FALSE, directed = FALSE, PoissonEdges = TRUE
  
  if(length(avgDeg)==0) return(fastRG(X,B, PoissonEdges=PoissonEdges, ...))  
  
  
  # if avgDeg is specified, then scale B by the appropriate amount.
  if(length(avgDeg) >0){
    eDbar = howManyEdges(X,B)[2]  # this returns the expected avg degree in the poisson graph.
    B = B * avgDeg/eDbar
  }
  
  
  if(!PoissonEdges){
    if(max(B)>=1){
      print(
        "This combination of B and avgDeg has led to probabilities that exceed 1.
        Suggestion:  Either diminish avgDeg or enable poisson edges."
        )
    }
      
    B = -log(1-B)   # this ensure that edge probabilites are bijand not 1-exp(-bij).
  }
  
  # avgDeg is set to NULL because it has already been handled internally.... this is because we need to scale B before making the transformation.
  return(fastRG(X,B, PoissonEdges = PoissonEdges, avgDeg=NULL,...))
}



fastRG <- function(X, S, Y= NULL, avgDeg = NULL,
                   simple = NULL, PoissonEdges = TRUE, directed = FALSE, selfLoops = FALSE){
  # X                    is an n x K1 matrix
  # S                    is a K1 x K2 matrix
  # Y                    is a d x K2 matrix; if Null, Y <- X
  # avgDeg               if specified and PoissonEdges == T, expected rowSums of output is avgDeg.
  #                      if specified and PoissonEdges == F, expected rowSums of output is less than avgDeg and close when output is sparse.
  # simple == T          sets PoissonEdges = FALSE, directed = FALSE, selfLoops = FALSE
  
  # selfLoops == T       retains the selfloops, 
  # directed == T        does not symmetrize the graph, 
  # PoissonEdges == T    keeps the multiple edges.
  
  # directed == F,       symmetrizes S and A.  avgDeg calculations are on out degree (i.e. rowSums)
  # selfLoops == F.      This uses a poisson approximation to the binomial...
  #                           Let M~poisson(\sum_{uv} \lambda_{uv}) be the number of edges.
  #                           If selfLoops == F, then the code uses the approximation 
  #                              Poisson(\lambda_{ij}) \\approx\\ Binomial(M, \lambda_{ij}/\sum_{uv}\lambda_{uv})
  #                              This approximation is good when total edges is ~n or larger and max \lambda_{ij} ~constant or smaller.
  # PoissonEdges == F    this thresholds each edge; multiple edges are replaced by single edges. 
  
  if(length(Y)>0){  # this means the output will be asymmetric and potentially rectangular
    directed = T
    selfLoops = T
    simple = F
  }
  if(is.null(Y)) Y <- X
  
  
  # Check that both X, Y, and S have non-negative entries.
  if(sum(X<0) + sum(S<0) + sum(Y<0)) return(NA)
  
  # if simple = T, then set PoissonEdges = FALSE, directed = FALSE, selfLoops = FALSE
  if(length(simple)>0) if(simple){
    selfLoops = FALSE
    directed= FALSE
    PoissonEdges = FALSE
  }
  
  
  n = nrow(X) 
  d = nrow(Y)
  K1 = ncol(X)
  K2 = ncol(Y)
  
  
  # if avgDeg is specified, then scale S by the appropriate amount.
  if(length(avgDeg) >0){
    eDbar = howManyEdges(X,S, Y)[2]  # this returns the expected avg degree in the poisson graph.
    S = S * avgDeg/eDbar
  }
  
  # if undirected, symmetrize S<- (S + t(S))/2 and then divide result by 2 because this doubles edge probabilities. 
  if(!directed) S = (S+t(S))/4
  
  
  
  Cx = diag(colSums(X))
  Cy = diag(colSums(Y))
  Xt  = apply(X,2,function(x) return(x/sum(x)))
  Yt  = apply(Y,2,function(x) return(x/sum(x)))
  St = Cx%*%S%*%Cy
  m = rpois(1, sum(St)) # number of edges to sample
  
  # if no edges, return empty matrix. 
  if (m==0) return(sparseMatrix(c(1:n),c(1:d),x = 0, dims = c(n, d)))
  
  # to sample U,V from St, we need to sample from a vector, then convert back to matrix.
  UV = sample(K1*K2, size = m, replace = TRUE, prob = St)  
  tabUV = table(c(1:(K1*K2), UV))
  Uindex = 1:(K1*K2) %% K1
  Uindex[Uindex==0]= K1
  Vindex = ((1:(K1*K2)) %/% K1 ) + 1
  Vindex = c(1, Vindex)[-(K1*K2 +1)]
  
  EdgeOut = rep(0, m)
  EdgeIn = EdgeOut
  ticker = 1
  
  for (uv in 1:(K1*K2)){
    edgesInThisBlock = tabUV[uv]
    I = sample(n, size = edgesInThisBlock, replace = T, prob = Xt[, Uindex[uv]])
    J = sample(d, size = edgesInThisBlock, replace = T, prob = Yt[, Vindex[uv]])
    
    EdgeOut[ticker:(ticker+edgesInThisBlock-1)] = I
    EdgeIn[ticker:(ticker+edgesInThisBlock-1)] = J
    ticker = ticker + edgesInThisBlock
  }
  
  A = sparseMatrix(EdgeOut, EdgeIn, x = 1, dims = c(n, d))
  if (selfLoops == FALSE){
    diag(A) = 0
  }
  if(!directed) A = A + t(A) # symmetrization doubles edge probabiilties! 
  if(!PoissonEdges) A@x[A@x>1]=1  # thresholding sets nonzero elements of A to one.
  
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




