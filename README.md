fastRG: a linear time algorithm for sampling a generalized random dot product graph.
==============================

fastRG quickly samples a generalized random dot product graph (RDPG), which is a generalization of a broad class of network models.   Given matrices $X \in R^{n \times K}$ and $S \in R^{K \times K}$ with positive entries, fastRG samples a matrix with  expectation $XSX^T$ and independent Poisson entries.  See Theorem \ref{theorem2} for an extension to Bernoulli entries.  See the forthcoming report for sampling Bernoulli entries. 

The basic idea of the algorithm is to first sample the number of edges, $m$, and then put down edges one-by-one.  It runs in $O(m)$ operations.  For example, in sparse graphs $m = O(n)$ and the algorithm is dramatically faster than ``element-wise'' algorithms which run in $O(n^2)$ operations.

Installation
------------

```R
source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
```



Functions 
------------
```R
fastRG(X, S, avgDeg = NULL, simple = NULL, PoissonEdges = TRUE, directed = FALSE, selfLoops = FALSE)
sbm(n,pi, B, PoissonEdges = F, ...)
dcsbm(theta,pi, B, ...)
dcMixed(theta,alpha, B, ...)
dcOverlapping(theta,pi, B, ...)

howManyEdges(X,S)
```
The functions sbm, dcsbm, dcMixed, and dcOverlapping are wrappers for fastRG.  

Arguments 
------------
```R
X              # X in the gRDPG
S              # S in the gRDPG
avgDeg         # to help ensure the graph is not too dense, avgDeg scales the S matrix to set 
#   the expected average expected degree to avgDeg. If avgDeg = null, this is ignored.

n              # number of nodes
pi             # a K vector of membership probabilities
alpha          # parameter of the dirichlet distribution in the assignment of block memberships in dcMixed
B              # middle probability matrix  (this becomes S in fastRG)
theta          # vector of degree parameter in degree corrected models, 
#  the expected adjacency matrix becomes: 
#            diag(theta) %*% X %*% S %*% t(X) %*% diag(theta)
#  in dcMixed and dcOverlapping, if theta is a single value, then 
#  n <- theta and there is no degree correction.


simple         # if TRUE, samples a simple graph by setting PoissonEdges = directed = multiEdges = FALSE
PoissonEdges   # See details
directed       # See details
selfLoops      # See details

```

Details
------------
fastRG samples a Bernoulli gRDPG where $\lambda_{ij} = X_i' S X_j$ and the probability of an edge $(i,j)$ is $1 - exp(-\lambda_{ij})$.  In sparse graphs, this is a good approximation to having edge probabilities $\lambda_{ij}$.  If multiEdges is set to TRUE, then it samples a Poisson gRDPG where $\lambda_{ij} = X_i' S X_j$.  Arugments can keep self loops or keep the graph directed.

sbm, dcsbm, dcOverlapping, and dcMixed are wrappers for fastRG that sample the Stochastic Blockmodel, Degree Corrected Stochastic Blockmodel, the Degree Corrected Overlapping Stochastic Blockmodel, and the Degree Corrected Mixed Membership Stochastic Blockmodel.  To remove Degree correction, set theta = rep(1, n).  

If selfLoops == T, then fastRG retains the selfloops. If selfLoops == F, then fastRG uses a poisson approximation to the binomial in the following sense: Let $M\sim poisson(\sum_{uv} \lambda_{uv})$ be the number of edges. fastRG approximates edge probabilities of $Poisson(\lambda_{ij})$ with  $Binomial(M, \lambda_{ij}/\sum_{uv}\lambda_{uv})$.  This approximation is good when total edges is order $n$ or larger and $\max \lambda_{ij}$ is order constant or smaller.


If directed == T, then fastRG does not symmetrize the graph.  If directed == F, then fastRG symmetrizes S and A.

If PoissonEdges == T, then fastRG keeps the multiple edges and avgDeg calculations are on out degree (i.e. rowSums).  If PoissonEdges == F, then fastRG thresholds each edge so that multiple edges are replaced by single edges. In this case, only SBM has edge probabilities exactly given by $\lambda$ (up to scaling for avgDeg).  The other techniques have edge probabilities $1 - exp(-\lambda)$.


Values
------------
fastRG and its wrappers all output a sparse matrix from the package Matrix. 

howManyEdges returns a vector with three elements.  The first element is the expected number of edges. The second is the expected average degree.  The third is the expected edge density. 

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

or fastRG also allows for simulating from E(A) = X S Y', where A and S could be rectangular.  This is helpful for bipartite graphs or matrices of features.

```R
n = 10000
d = 1000
K1 = 5
K2 = 3

X = matrix(rpois(n = n*K1, 1), nrow = n)
Y = matrix(rpois(n = d*K2, 1), nrow = d)
S = matrix(runif(n = K1*K2, 0,.1), nrow = K1)

A = fastRG(X,S,Y, avgDeg = 10)
```





or

```R
K = 10
n = 500
pi = rexp(K) +1
pi = pi/sum(pi) * 3
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
theta = rexp(n)

A = dcsbm(theta, pi,B,avgDeg = 50)
# here is the average degree:
mean(rowSums(A))  

# If we remove multiple edges, the avgDeg parameter is not trustworthy:
A = dcsbm(theta, pi,B,avgDeg = 50, PoissonEdges = F)
mean(rowSums(A))  

# but it is a good upper bound when the graph is sparse:
n = 10000
A = dcsbm(rexp(n), pi,B,avgDeg = 50, PoissonEdges = F)
mean(rowSums(A))  

```      

or

```R
# This draws a 500 x 500 adjacency matrix from each model.
#   Each image might take around 5 seconds to render.

K = 10
n = 500
pi = rexp(K) +1
pi = pi/sum(pi) * 3
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
theta = rexp(n)
A= dcsbm(theta, pi,B,avgDeg = 50)
image(as.matrix(t(A[,n:1])),col = grey(seq(1,0, len=20)))


K = 2
n = 500
alpha = c(1,1)/5
B = diag(c(1,1))
theta = n
A= dcMixed(theta, alpha,B,avgDeg = 50)
image(as.matrix(t(A[,theta:1]))/max(A),col = grey(seq(1,0, len=20)))


n = 500
K = 2
pi = c(.7,.7)
B = diag(c(1,1))
theta = n
A= dcOverlapping(theta, pi,B,avgDeg = 50)
image(as.matrix(t(A[,n:1]))/max(A),col = grey(seq(1,0, len=20)))


K = 10
n = 500
pi = rexp(K) +1
pi = pi/sum(pi) 
B = matrix(rexp(K^2), nrow=K) 
B = B/ (3*max(B))
diag(B) = diag(B)+ mean(B)*3
A= sbm(n, pi,B)
image(as.matrix(t(A[,n:1])),col = grey(seq(1,0, len=20)))
mean(A)
```




```R
# this samples a DC-SBM with 10,000 nodes
#  then computes, and plots the leading eigenspace.
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
A = dcsbm(rgamma(n,shape = 2,scale = .4), pi,B,avgDeg = 20, simple = T)
mean(rowSums(A))

# leading eigen of regularized Laplacian with tau = 1
D = Diagonal(n, 1/sqrt(rowSums(A)+1))
ei = eigs_sym(D%*%A%*%D, 10)  

# normalize the rows of X:
X = t(apply(ei$vec[,1:K],1, function(x) return(x/sqrt(sum(x^2)+1/n))))

# taking a varimax rotation makes the leading vectors pick out clusters:
X = varimax(X, normalize = F)$load

par(mfrow = c(5,1), mar = c(1,2,2,2), xaxt = "n",yaxt = "n")

# plot 1000 elements of the leading eigenvectors:
s = sort(sample(n,1000))
for(i in 1:5){
plot(X[s,i], pch  ='.')
}
dev.off()

```    

or



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
ei = eigs_sym(L, 4)  

s = sort(sample(n, 10000))
X = t(apply(ei$vec[,1:K],1, function(x) return(x/sqrt(sum(x^2)+1/n))))
plot(X[s,3])  # highly localized eigenvectors
```   



To sample from a degree corrected and node contextualized graph...

```R
n = 10000  # number of nodes
d = 1000 # number of features
K = 5  # number of blocks 


# to generate an X for dc-sbm, we need pi and theta.

pi = rexp(K) +1
pi = pi/sum(pi) * 3
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
theta = rep(0,n)

# to make pretty pictures, order by cluster size:
B = B[order(pi), ]
B = B[, order(pi)]
pi = sort(pi/sum(pi))

z = sample(K,n,replace = T, prob = pi)
# Again, so that it makes pictures...
z = sort(z)  
X = sparse.model.matrix(~as.factor(z)-1)
Tz = table(z)
ct= c(0,cumsum(Tz))
for(i in 1:K){
  theta[(ct[i]+1):ct[i+1]] = -sort(-rexp(Tz[i]))
}
X@x = theta

# first the graph:
A = fastRG(X,B, avgDeg = 10)

# now the features:
X@x = X@x + rexp(n) # the degree parameter should be different. X@x + rexp(n) makes feature degrees correlated to degrees in graph.


thetaY = rep(0, d)
piFeatures = rexp(K) +1
piFeatures = piFeatures/sum(piFeatures) * 3
BFeatures = matrix(rexp(K^2)+1, nrow=K)
diag(BFeatures) = diag(BFeatures)+ mean(BFeatures)*K

# to make pretty pictures, order by cluster size:
BFeatures = BFeatures[, order(piFeatures)]
piFeatures = sort(piFeatures/sum(piFeatures))

y = sample(K,d,replace = T, prob = piFeatures)
# Again, so that it makes pictures...
y = sort(y)  
Y = sparse.model.matrix(~as.factor(y)-1)
Ty = table(y)
cty= c(0,cumsum(Ty))
for(i in 1:K){
  thetaY[(cty[i]+1):cty[i+1]] = -sort(-rexp(Ty[i]))
}
Y@x = thetaY

# now generate the features:

features = fastRG(X,BFeatures, Y, avgDeg = 20)

```


