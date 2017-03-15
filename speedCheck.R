# this code generates a few 5M edge graphs and compares various processing times.

source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
n = 500000
K = 3

X = matrix(rpois(n = n*K, 1), nrow = n)
S = matrix(runif(n = K*K, 0,.0001), nrow = K)

# when sampling a fastRG, the generation of the edge list is often faster than
#   the time it takes to convert the edge list into an igraph:

library(igraph)
# first, just make the edgelist:
system.time({el = fastRG(X,S,avgDeg = 10, returnEdgeList = T)})

# form an igraph from an edge list:
#  (compare this speed to the last line)
system.time({g = graph_from_edgelist(el)})


# forming the edge list is also faster than the time it takes to turn the edge list into a sparse Matrix:
# first, just make the edgelist:
system.time({el = fastRG(X,S,avgDeg = 10, returnEdgeList = T)})

# this next line forms the edge list then converts to a sparse Matrix.
#  (this takes more than twice as long):
system.time({A = fastRG(X,S, avgDeg = 10)})

# it is faster to create the pattern Matrix (i.e. ignore repeated edges multiple):
system.time({A = fastRG(X,S, PoissonEdges = F, avgDeg = 10)})



#  it is possible that different parameterizations of the fastRG lead to different results.