fastRG: a linear time algorithm for sampling a generalized random product graph.
==============================

fastRG quickly samples a generalized random product graph, which is a generalization of a broad class of network models.   Given matrices $X \in R^{n \times K_x}$, $Y \in R^{n \times K_y}$, and $S \in R^{K_x \times K_y}$ each with positive entries, fastRG samples a matrix with  expectation $XSY^T$ and independent Poisson entries.  [See the tech report for sampling Bernoulli entries](https://arxiv.org/abs/1703.02998). 

The basic idea of the algorithm is to first sample the number of edges, $m$, and then put down edges one-by-one.  It runs in $O(m)$ operations.  For example, in sparse graphs $m = O(n)$ and the algorithm is dramatically faster than ``element-wise'' algorithms which run in $O(n^2)$ operations.

Installation in R
------------

```R
source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
```



Functions 
------------
```R
fastRG(X, S, Y= NULL, avgDeg = NULL, simple = NULL, 
          PoissonEdges = TRUE, directed = FALSE, selfLoops = FALSE, 
          returnEdgeList = FALSE, returnParameters = FALSE)
er(n, p = NULL, avgDeg =NULL, directed = FALSE, returnEdgeList = FALSE,...)     
cl = function(theta, avgDeg = NULL, directed = FALSE, returnEdgeList = FALSE,...)   
sbm(n,pi, B, PoissonEdges = F, returnParameters = FALSE, parametersOnly = FALSE, ...)
dcsbm(theta,pi, B, returnParameters = FALSE, parametersOnly = FALSE, ...)
dcMixed(theta,alpha, B, returnParameters = FALSE, parametersOnly = FALSE, ...)
dcOverlapping(theta,pi, B, returnParameters = FALSE, parametersOnly = FALSE, ...)

howManyEdges(X,S)
```
The functions er, sbm, dcsbm, dcMixed, and dcOverlapping are wrappers for fastRG. 

Arguments 
------------
```R
X                 # X in the gRDPG
S                 # S in the gRDPG

Y                 # if Null, Y <- X; if not, E(A) = XSY', directed <- T, selfLoops <- T, simple <- F
                  #   Y need not have same number of rows or columns as X.  matrix mult X %*% S %*% t(Y) must be defined.
  
avgDeg            # to help ensure the graph is not too dense, avgDeg scales the S matrix to set 
                  #   the expected average expected degree to avgDeg. If avgDeg = null, this is ignored.

n                 # number of nodes
pi                # a K vector of membership probabilities
alpha             # parameter of the dirichlet distribution in the assignment of block memberships in dcMixed
B                 # middle probability matrix  (this becomes S in fastRG)
theta             # vector of degree parameter in degree corrected and Chung-Lu models, 
                  #  the expected adjacency matrix becomes: 
                  #           diag(theta) %*% X %*% S %*% t(X) %*% diag(theta)
                  #  in dcMixed and dcOverlapping, if theta is a single value, then 
                  #  n <- theta and there is no degree correction.
returnParameters  # if TRUE, then it returns a list(A, X, S, Y).  This can be helpful for simulation studies which 
                              # seek to estimate X, S, Y.  If FALSE, then it returns A.  
                              # if returnEdgeList = TRUE, then this parameter is ignored.
parametersOnly    # if TRUE, then the wrapper only returns the X and S matrix that would otherwise be sent to fastRG.

simple            # if TRUE, samples a simple graph by setting PoissonEdges = directed = multiEdges = FALSE
PoissonEdges      # See Details
directed          # See Details
selfLoops         # See Details
returnEdgeList    # See Values
```

Details
------------
fastRG samples a Poisson gRPG where $\lambda_{ij} = X_i' S Y_j$ is the rate parameter for edge $i,j$.   If multiEdges is set to FALSE, then it samples a Bernoulli gRPG where the probability of edge $(i,j)$ is $1 - exp(-\lambda_{ij})$.  In sparse graphs, this is a good approximation to having edge probabilities $\lambda_{ij}$. Arugments can keep self loops or keep the graph directed.

er, cl, sbm, dcsbm, dcOverlapping, and dcMixed are wrappers for fastRG that sample the Erdos-Renyi, Chung-Lu, Stochastic Blockmodel, Degree Corrected Stochastic Blockmodel, the Degree Corrected Overlapping Stochastic Blockmodel, and the Degree Corrected Mixed Membership Stochastic Blockmodel.  To remove Degree correction, set theta = rep(1, n) or set theta equal to the number of desired nodes $n$.

If selfLoops == T, then fastRG retains the selfloops. If selfLoops == F, then fastRG uses a poisson approximation to the binomial in the following sense: Let $M\sim poisson(\sum_{uv} \lambda_{uv})$ be the number of edges. fastRG approximates edge probabilities of $Poisson(\lambda_{ij})$ with  $Binomial(M, \lambda_{ij}/\sum_{uv}\lambda_{uv})$.  This approximation is good when total edges is order $n$ or larger and $\max \lambda_{ij}$ is order constant or smaller.  Under er and sbm, there is a default correction that removes this issue. 

If directed == T, then fastRG does not symmetrize the graph.  If directed == F, then fastRG symmetrizes S and A.

If PoissonEdges == T, then fastRG keeps the multiple edges and avgDeg calculations are on out degree (i.e. rowSums).  If PoissonEdges == F, then fastRG thresholds each edge so that multiple edges are replaced by single edges. In this case, only SBM has edge probabilities exactly given by $\lambda$ (up to scaling for avgDeg).  The other techniques have edge probabilities $1 - exp(-\lambda)$.

If Y is specified, then it returns a sparse matrix, with poisson entries, where $E(A) = X S Y'$.  This is useful for creating sparse rectangular matrices.  For example, this could be a feature matrix.  Or, it could be a "rectangular adjacency matrix" for a bipartite graph.


Values
------------
By default, fastRG and its wrappers all output a sparse matrix from the package Matrix. If returnEdgeList = TRUE, then it returns an edge list matrix with m rows and 2 columns. 

howManyEdges returns a vector with two elements.  The first element is the expected number of edges. The second is the expected average degree.  

Example Usage
-------------
First, generate the edgelist for an Erdos-Renyi graph with n = 1,000,000 nodes and expected degree 5.  That makes the edge probability 5/n.  
```R
# install fastRG:
source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
# sample:
system.time(er(n=10^6, avgDeg = 5, returnEdgeList = F))
```

```R
n = 10000
K = 5

X = matrix(rpois(n = n*K, 1), nrow = n)
S = matrix(runif(n = K*K, 0,.0001), nrow = K)

howManyEdges(X,S)[-1]
A = fastRG(X,S, simple=T)

# if you want to create an igraph, it is fastest to
#   1) return an edgelist from fastRG and 
#   2) form the igraph from the edgelist.

library(igraph)
el = fastRG(X,S, simple=T, returnEdgeList = T)
g = graph_from_edgelist(el)
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
# This draws a 100 x 100 adjacency matrix from each model.
#   Each image might take around 5 seconds to render.

K = 10
n = 100
pi = rexp(K) +1
pi = pi/sum(pi) * 3
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
theta = rexp(n)
A= dcsbm(theta, pi,B,avgDeg = 50)
image(as.matrix(t(A[,n:1])),col = grey(seq(1,0, len=20)))


K = 2
n = 100
alpha = c(1,1)/5
B = diag(c(1,1))
theta = n
A= dcMixed(theta, alpha,B,avgDeg = 50)
image(as.matrix(t(A[,theta:1]))/max(A),col = grey(seq(1,0, len=20)))


n = 100
K = 2
pi = c(.7,.7)
B = diag(c(1,1))
theta = n
A= dcOverlapping(theta, pi,B,avgDeg = 50)
image(as.matrix(t(A[,n:1]))/max(A),col = grey(seq(1,0, len=20)))


K = 10
n = 100
pi = rexp(K) +1
pi = pi/sum(pi) 
B = matrix(rexp(K^2), nrow=K) 
B = B/ (3*max(B))
diag(B) = diag(B)+ mean(B)*K
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

A <- fastRG::dcsbm(theta = rgamma(n,shape = 2,scale = .4), pi = pi, B = B, avg_deg = 10)

A = dcsbm(rgamma(n,shape = 2,scale = .4), pi,B,avgDeg = 10)

set.seed(39)
theta <- rgamma(n,shape = 2,scale = .4)

set.seed(40)
p = dcsbm(theta, pi,B,avgDeg = 10, parametersOnly = T)

set.seed(40)
p2 = dcsbm_params(theta, pi,B, avg_deg = 10)

all.equal(p$X, p2$X)
all.equal(p$S, p2$S)

e <- howManyEdges(p$X, p$S)
e2 <- expected(p$X, p$S)

all.equal(e[1], e2[[1]])
all.equal(e[2], e2[[2]])

fastRG::fastRG(p$X, p$S)
fastRG(p$X, p$S)




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


# Here are the parameters for the graph:

pi = rexp(K) +1
pi = pi/sum(pi) * 3
B = matrix(rexp(K^2)+1, nrow=K)
diag(B) = diag(B)+ mean(B)*K
theta = rexp(n)
paraG = dcsbm(theta=theta, pi = pi, B=B,parametersOnly = T)


# Here are the parameters for the features:

thetaY = rexp(d)
piFeatures = rexp(K) +1
piFeatures = piFeatures/sum(piFeatures) * 3
BFeatures = matrix(rexp(K^2)+1, nrow=K)
diag(BFeatures) = diag(BFeatures)+ mean(BFeatures)*K

paraFeat = dcsbm(theta = thetaY,pi = piFeatures, B = BFeatures,parametersOnly = T)

# the node "degrees" in the features, should be related to their degrees in the graph.
X = paraG$X
X@x = paraG$X@x + rexp(n) # the degree parameter should be different. X@x + rexp(n) makes feature degrees correlated to degrees in graph.



# generate the graph and features
A = fastRG(paraG$X,paraG$S, avgDeg = 10)
features = fastRG(X,paraFeat$S, paraFeat$X, avgDeg = 20)
```



This next bit of code repeats the computational experiment from Section 4.1.


```R
library(tidyverse)

logSeq = function(from, to, len){
  # find a sequence from, ..., to of length len such that *on the log scale* the points are equidistant. 
  seq(log(from), log(to), len =len) %>% exp %>%  round %>% return
}


nseq = logSeq(10000, 10000000,len = 10)
mseq =  logSeq(100000, 100000000,len = 10)


K = 5
runTimes = matrix(NA, nrow = length(mseq)*length(nseq), ncol=3) %>% as_data_frame %>% as.tbl
colnames(runTimes) = c("n","m","time")

S = matrix(runif(n = K*K), nrow = K)
ticker = 1
for(ntick in length(nseq):1){
  n = nseq[ntick]
  X = matrix(rpois(n = n*K, 1), nrow = n)
  
  for(mtick in 1:length(mseq)){
    averageDegree = mseq[mtick]/n
    rm(el)
    gc()
    
    timer = system.time({
      el = fastRG(X,S, avgDeg = averageDegree, PoissonEdges = T,directed = T,selfLoops = T, returnEdgeList= T)
    })
    
    runTimes[ticker,] = rbind(n,mseq[mtick],timer[3])
    ticker  = ticker+1
  }
}




pdf(file = "runTimeFixN.pdf", height =4, width = 4)
runTimes %>% mutate("log10(n)" = as.character(round(log10(n),2)), "E(m)"= m) %>% 
  ggplot(aes(x=`E(m)`, y=time, group = `log10(n)`))+
  geom_line(aes(color = `log10(n)`)) + guides(col = guide_legend(reverse = TRUE))+
  geom_point()+
  scale_x_log10() +scale_y_log10() + 
  geom_abline(slope  = 1, intercept = -6.85,  size = 1) + 
  ggtitle("Fixing n, the running time appears\nto grow linearly with E(m)")
dev.off()

pdf(file = "runTimeFixM.pdf", height =4, width = 4)
runTimes %>% mutate("log10(E(m))" = as.character(round(log10(m),2))) %>% 
  ggplot(aes(x=n, y=time, group = `log10(E(m))`))+
  geom_line(aes(color = `log10(E(m))`)) + guides(col = guide_legend(reverse = TRUE))+
  geom_point()+
  scale_x_log10() +scale_y_log10() + 
  geom_abline(slope  = 1, intercept = -6.3, size =1)+
  ggtitle("Fixing E(m), the running time appears\nto grow linearly with n")
dev.off()



```
