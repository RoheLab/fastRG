# Fast way of sampling Random dot product graph
# Karl Rohe, Jun Tao, Xintian Han, Norbert Binkiewicz 2017
# For details of the functions see https://github.com/karlrohe/fastRG/blob/master/README.md
# based upon this paper:  https://arxiv.org/abs/1703.02998

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
  # density = em/(nrow(X)*nrow(Y))  # in big graphs, nrow(X)*nrow(Y) causes integer overflow errors...
  return(c(em,avDeg))
  # return(c(em,avDeg, density))
  
}


dcOverlapping = function(theta,pi, B, returnParameters = FALSE, parametersOnly = FALSE, ...){
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
  
  if(parametersOnly) return(list(X=X, S=B))
  return(fastRG(X,B, returnParameters = returnParameters, ...))
}



dcMixed = function(theta,alpha, B, returnParameters = FALSE, parametersOnly = FALSE, ...){
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
  
  if(parametersOnly) return(list(X=X, S=B))
  return(fastRG(X,B, returnParameters = returnParameters, ...))
}



dcsbm = function(theta,pi, B, returnParameters = FALSE, parametersOnly = FALSE, ...){
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
  
  if(parametersOnly) return(list(X=X, S=B))
  return(fastRG(X,B, returnParameters = returnParameters, ...))
  
}










sbm = function(n,pi, B, PoissonEdges = F, avgDeg =NULL, returnParameters = FALSE, parametersOnly = FALSE, ...){
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
  
  if(length(avgDeg)==0){
    return(fastRG(X,B, PoissonEdges=PoissonEdges, ...))  
  }
  
  
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
  if(parametersOnly) return(list(X=X, S=B))
  return(fastRG(X,B, PoissonEdges = PoissonEdges, avgDeg=NULL, returnParameters = returnParameters, ...))
}



fastRG <- function(X, S, Y= NULL, avgDeg = NULL,
                   simple = NULL, PoissonEdges = TRUE, directed = FALSE, 
                   selfLoops = FALSE, returnEdgeList = FALSE, returnParameters = FALSE){
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
    returnY = TRUE
  }
  if(is.null(Y)){
    Y <- X
    returnY = FALSE
  }
  
  
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
  # Xt  = apply(X,2,function(x) return(x/sum(x)))
  # Yt  = apply(Y,2,function(x) return(x/sum(x)))  # sample(n,...) automatically does the normalization.
  St = Cx%*%S%*%Cy
  m = rpois(n = 1, lambda = sum(St)) # number of edges to sample
  
  # if no edges, return empty matrix. 
  if (m==0){
    A = sparseMatrix(c(1:n),c(1:d),x = 0, dims = c(n, d))
    if(returnParameters){
      if(returnY) out = list(A = A, X = X, S = S, Y = Y)
      if(!returnY) out = list(A = A, X = X, S = S, Y = NULL)
    }
    if(!returnParameters){
      out = A
    }
    return(out)
  }
  
  # this simulates \varpi, denoted here as tabUV.  element u,v is the number of edges between column u and column v.
  tabUV = matrix(rmultinom(n=1, size= m, prob = St), nrow = K1, ncol = K2)
  cumsumUV = matrix(cumsum(tabUV), nrow = K1, ncol = K2)
  
  #  cbind(eo,ei) is going to be the edge list.  eo = "edge out node".  ei = "edge in node"
  eo = rep(NA, m)  
  ei = eo
  
  
  # to avoid doing K1*K2 samples for I and J, we can instead take
  #   only K1 + K2 samples.  This requires some awkward indexing.
  #   for(u in 1:K1) for(v in 1:K2) tabUV[u,v]  <- eventually this is going to happen.
  #     because this first sets u =1 and loops through the different values of v, 
  #     the awkward indexing will only apply to the "v" or the ei vector.
  #     for eo, we can just loop through u in 1:K1...
  
  
  ticker = 1
  blockDegreesU = rowSums(tabUV)
  for(u in 1:K1){
    if(blockDegreesU[u]>0){
      eo[ticker:(ticker + blockDegreesU[u]-1)] = 
        sample(n, size = blockDegreesU[u], replace = T, prob = X[, u])
      ticker = ticker + blockDegreesU[u]
    }
  }
  
  
  # for ei, things are more awkward.  if we had pointers, perhaps there would be a faster way... instead
  #   create ei-tmp.  which will hold the values of ei, but in the wrong order.  another loop will take care of the indexing.
  eitmp = ei  
  ticker = 1
  blockDegreesV = colSums(tabUV)
  for(v in 1:K2){
    if(blockDegreesV[v]>0){
      eitmp[ticker:(ticker + blockDegreesV[v]-1)] = sample(d, size = blockDegreesV[v], replace = T, prob = Y[, v])
      ticker = ticker + blockDegreesV[v]
    }
  }
  
  
  #  this loop worries about the indexing... putting ei-tmp in the correct order, to match up with eo. 
  ticker = 1
  # tickerU = 1
  tickerV = c(1,cumsum(blockDegreesV))
  for (u in 1:K1){ 
    for(v in 1:K2){
      edgesInThisBlock = tabUV[u,v]
      if(edgesInThisBlock>0){
        # I = sample(n, size = edgesInThisBlock, replace = T, prob = X[, Uindex[uv]])
        # J = sample(d, size = edgesInThisBlock, replace = T, prob = Y[, Vindex[uv]])
        
        # eo[ticker:(ticker+edgesInThisBlock-1)] = eotmp[tickerU:(tickerU+edgesInThisBlock-1)]
        # tickerU = tickerU + edgesInThisBlock
        ei[ticker:(ticker+edgesInThisBlock-1)] = eitmp[tickerV[v]:(tickerV[v]+edgesInThisBlock-1)]
        tickerV[v] = tickerV[v] + edgesInThisBlock
        ticker = ticker + edgesInThisBlock
      }
    }
  }
  
  
  if (!selfLoops){ 
    goodEdges = eo!=ei
    eo = eo[goodEdges]
    ei =  ei[goodEdges]
  }
  
  if(!directed){  if(n!=d) break
    # if it is directed and has self-loops, the current implementation has the corrected expected number of edges... but it is required to be an even number... not poisson.
    # symmetrization doubles edge probabiilties! 
    eoOLD = eo
    eo = c(eo, ei)
    ei = c(ei, eoOLD)
  }
  
  if(returnEdgeList) return(cbind(eo,ei))
  
  if(PoissonEdges) A = sparseMatrix(eo, ei, x = 1, dims = c(n, d))
  if(!PoissonEdges) A = sparseMatrix(i = eo, j = ei,dims = c(n, d))  # thresholding sets nonzero elements of A to one.
  
  if(returnParameters){
    if(returnY) out = list(A = A, X = X, S = S, Y = Y)
    if(!returnY) out = list(A = A, X = X, S = S, Y = NULL)
    return(out)
  }
  if(!returnParameters){
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




