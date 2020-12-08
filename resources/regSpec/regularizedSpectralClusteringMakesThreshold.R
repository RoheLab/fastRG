# how does reg spec perform at reconstruction threshold?


source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
library(rARPACK)

n  =10000  # number of nodes in graph.
replicates = 20  # how many times to simulate each graph?


# here is a single simulation.
a = 4.5  # this is expected "in block degree"
b = 1  # this is expected "out block degree"
s = (a - b)/2
p = (a + b)/2
# if this is greater than 1, then we are above reconstruction threshold.
s^2 / p

# sample a graph
B = matrix(c(a,b,b,a), nrow=2)/n
A = sbm(n,c(1,1), B)

# perform regularized spectral clustering with tau = mean degree.
rs = rowSums(A)
D = Diagonal(n, 1/sqrt(rs + mean(rs)))
L = D%*%A%*%D
ei = eigs(L, 10)
v = ei$vectors[,2]

# here is a scree plot:
plot(abs(ei$values[-1]))
# here is the eigenvector
plot(v, pch = ".")

# above reconstruction threshold, it should have positive correlation with true partition:
abs(cor(v, c(rep(1,n/2), rep(0,n/2))))


# now lets do that a bunch of times.


# this function performs the simulation for a,b,n.
sim = function(a,b,n){
  # given a,b, and n.
  # simulates SBM.  computes 2nd eigen of reg laplacian and returns its correlation with true.
  B = matrix(c(a,b,b,a), nrow=2)/n
  A = sbm(n,c(1,1), B)
  rs = rowSums(A)
  D = Diagonal(n, 1/sqrt(rs + mean(rs)))
  L = D%*%A%*%D
  v = eigs(L, 2)$vec[,2]
  # return(cor.test(v, c(rep(1,n/2), rep(0,n/2)),method = "kend")$p.value)
  return(abs(cor(v, c(rep(1,n/2), rep(0,n/2)))))
}


# now, for a range of a,b values, simulate "replicates" many graphs.
# to change n, go back up to the top.  or set it here
#  n= 10000
#  replicates = 10
cc = matrix(NA, nrow = 30, ncol = replicates)
for(btick in 1:30){
  a = 5
  b = seq(1,1.75, len = 30)[btick]
  for(rep in 1:ncol(cc)) cc[btick,rep] = sim(a,b,n)
}

# now, plot the results.
a = 5
b = seq(1,1.75, len = 30)
s = (a - b)/2
p = (a + b)/2
# s^2; p
pdf(file ="regSpec/Reconstruction.pdf")
plot(s^2/p, rowMeans(cc), main= "when s^2/p >1, then the blocks can be reconstructed.\nif reg spec works, then eigenvector should correlate in that regime", ylab= "correlation of eigenvector with true partition")
lines(c(1,1), c(-10,10))
dev.off()

