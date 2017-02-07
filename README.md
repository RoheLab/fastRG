fastRG: a linear time algorithm for sampling a generalized random dot product graph.
==============================

fastRG quickly samples a generalized random dot product graph (RDPG), which is a generalization of a broad class of network models.   Given matrices $X \in R^{n \times K}$ and $S \in R^{K \times K}$ with positive entries, fastRG samples a matrix with  expectation $XSX^T$ and independent Poisson entries.  See Theorem \ref{theorem2} for an extension to Bernoulli entries.  See the forthcoming report for sampling Bernoulli entries. 

The basic idea of the algorithm is to first sample the number of edges, $m$, and then put down edges one-by-one.  It runs in $O(m)$ operations.  For example, in sparse graphs $m = O(n)$ and the algorithm is dramatically faster than ``element-wise'' algorithms which run in $O(n^2)$ operations.

Given $X$ and $S$, howManyEdges(X,S) returns the expected number of edges, the expected average node degree, and the expected edge density.  

Usage
------------
```R
sbm(n,pi, B, avgDeg=10)
dcsbm(theta,pi, B, avgDeg=10)
dcOverlapping(theta,pi, B, avgDeg=10)
dcMixed(theta,alpha, B, avgDeg=10)
fastRG(X, S,simple = NULL, selfLoops = FALSE, directed = FALSE, multiEdges = TRUE)
howManyEdges(X,S)
```

Arguments 
------------
```R
n          # number of nodes
pi         # a K vector of membership probabilities
B          # middle probability matrix
theta      # degree parameter in degree corrected models
alpha      # parameter of the dirichlet distribution in the assignment of block memberships.
avgDeg     # to help ensure the graph is not too dense, avgDeg scales the B matrix to set 
           #   the expected average degre. If avgDeg = null, this is ignored.
X          # X in the gRDPG
S          # S in the gRDPG
simple     # if TRUE, then sets selfLoops = directed = multiEdges = FALSE
```

Details
------------
fastRG samples a Bernoulli gRDPG where $\lambda_{ij} = X_i' S X_j$ and the probability of an edge $(i,j)$ is $1 - exp(-\lambda_{ij})$.  In sparse graphs, this is a good approximation to having edge probabilities $\lambda_{ij}$.  If multiEdges is set to TRUE, then it samples a Poisson gRDPG where $\lambda_{ij} = X_i' S X_j$.  Arugments can keep self loops or keep the graph directed.

sbm, dcsbm, dcOverlapping, and dcMixed are wrappers for fastRG that sample the Stochastic Blockmodel, Degree Corrected Stochastic Blockmodel, the Degree Corrected Overlapping Stochastic Blockmodel, and the Degree Corrected Mixed Membership Stochastic Blockmodel.  To remove Degree correction, set theta = rep(1, n).  

howManyEdges returns a vector with three elements.  The first element is the expected number of edges. The second is the expected average degree.  The third is the expected edge density. 

Note:  Only SBM has edge probabilities exactly given by $\lambda$ (up to scaling for avgDeg).  The other techniques have edge probabilities $1 - exp(-\lambda)$


Installation
------------

```R
source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
```

Example Usage
-------------

```R
n = 10000
K = 5

X = matrix(rpois(n = n*K, 1), nrow = n)
S = matrix(runif(n = K*K, 0,.0001), nrow = K)

howManyEdges(X,S)[-1]
A = fastRG(X,S, simple=T)
```

or

```R
K = 10
n = 500
pi = rexp(K) +1
pi = pi/sum(pi) * 3
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
A = dcsbm(rgamma(n,shape = 2,scale = .4), pi,B,avgDeg = 20)
hist(rowSums(A))
# average degree is a bit smaller than 10 because we have simplified the graph, removing multiple edges:
mean(rowSums(A))  
image(as.matrix(t(A[,n:1])),col = grey(c(1,0)))
```


```R
# this samples a DC-SBM with 10,000 nodes
#  then finds, and plots the leading eigenspace.
#  the code should run in less than a second.

require(rARPACK)
require(stats)
K = 10
n = 10000
pi = rexp(K) +1
pi = pi/sum(pi) 
pi = -sort(-pi)
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
A = dcsbm(rgamma(n,shape = 2,scale = .4), pi,B,avgDeg = 20)

# leading eigen of regularized Laplacian with tau = 1
D = Diagonal(n, 1/sqrt(rowSums(A)+1))
ei = eigs_sym(D%*%A%*%D, 10)  
# normalize the rows of X:
X = t(apply(ei$vec[,1:K],1, function(x) return(x/sqrt(sum(x^2)+1/n))))
# taking a varimax rotation makes the leading vectors pick out clusters:
X = varimax(X, normalize = F)$load
par(mfrow = c(5,1), mar = c(0,1,1,0), bty = "n", xaxt = "n")
# plot leading eigenvectors:
for(i in 1:5){
  plot(X[,i], pch  ='.')
  lines(c(0,n), c(0,0))
}
dev.off()

```



```R
# This samples a 1M node graph.  
# Depending on the computer, sampling the graph should take between 10 and 30 seconds 
#   Then, taking the eigendecomposition of the regularized graph laplacian should take between 1 and 3 minutes
# The resulting adjacency matrix is a bit larger than 100MB.
# The leading eigenvectors of A are highly localized
K = 3
n = 1000000
pi = rexp(K) +1
pi = pi/sum(pi) 
pi = -sort(-pi)
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
A = dcsbm(rgamma(n,shape = 2,scale = .4), pi,B,avgDeg = 10)
D = Diagonal(n, 1/sqrt(rowSums(A)+10))
L = D%*%A%*%D
ei = eigs_sym(L, 10)  
s = sort(sample(n, 10000))
X = t(apply(ei$vec[,1:K],1, function(x) return(x/sqrt(sum(x^2)+1/n))))
plot(X[s,2])  # highly localized eigenvectors

Ac=A
core = 3
while(min(rs) <= core){
good = rs>core
  Ac = Ac[good,good]
  rs = rowSums(Ac)
}
dim(Ac)
ei = eigs_sym(Ac, 4)  
s = sort(sample(nrow(Ac), 10000))
plot(ei$vec[s,2])  # highly localized eigenvectors


```