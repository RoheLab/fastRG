#' Sample a random dot product graph (RDPG)
#'
#' @param X An `n` by `k_1` matrix.
#' @param S A `k_1` by `k_2` matrix.
#' @param Y A `d` by `k_2` matrix. Defaults to `X`.
#' @param avgDeg When specified, the expected average degree of nodes
#'   in the output poisson multi-graph. When `poissonEdges = FALSE`, the
#'   resulting graph will have lower average degree due to lack
#'   of multiple-edges. When the graph is sparse, the expected number of
#'   edges for the Poisson multi-graph and Bernoulli graph are nearly the
#'   same. Defaults to `NULL`, such that no scaling occurs.
#' @param simple When `TRUE` indicates that you want to work with undirected
#'   graphs where self-loops and multi-edges are prohibited. Accomplishes
#'   this by setting `directed = FALSE`, `selfLoops = FALSE`, and
#'   `PoissonEdges = FALSE`, and then ignoring arguments `directed`,
#'   `selfLoops` and `PoissonEdges`. Defaults to `FALSE`.
#' @param PoissonEdges Logical indicating whether or not multi-edges are
#'   allowed. Defaults to `TRUE`, which keeps multi-edges and produces
#'   a multi-graph. When `FALSE`, only single edges are allowed, resulting
#'   in a graph. See **details** for some additional comments. Effected by
#'   `simple` argument.
#' @param directed Logical indicating whether or not the graph should be
#'   directed. Default is `directed = TRUE`. When `directed = FALSE`,
#'   symmetrizes `S` internally. When `X = Y` (which is the default when
#'   no `Y` is specified), this results in a symmetric adjacency matrix
#'   as output. When `avgDeg` is specified and the desired graph is directed,
#'   the average degree scaling is on the out-degree of each node (or the
#'   row sums if you prefer to think in terms of the adjacency matrix).
#'   Effected by the `simple` argument.
#' @param selfLoops Logical indicating whether edges are allowed from
#'   a node back to itself. Defaults to `TRUE`. When `FALSE`, sampling
#'   proceeds normally, and then self-loops are removed after sampling
#'   is completed. Effected by the `simple` argument.
#' @param returnEdgeList Logical indicating whether to return an edgelist
#'   rather than an adjacency matrix.
#'
#' @return A random Poisson (or Bernoulli) dot product graph. By default,
#'   returns a [Matrix::sparseMatrix()] in CSC format (i.e. of abstract
#'   class `CsparseMatrix`). When the graph is undirected, the `sparseMatrix`
#'   will also be symmetric.
#'
#'   If `returnEdgeList = TRUE`, instead returns the graph as an edgelist
#'   in matrix format, where the first column of the matrix contains two
#'   columns `from` and `to` . The `from` column contains the index
#'   of the outgoing edge in this case, and the `to` column contains
#'   the index of the incoming edge.
#'
#' @details TODO: clean up the following:
#'
#'   This uses a poisson approximation to the binomial...
#'
#'   Let \eqn{M ~ Poisson(\sum_{uv} \lambda_{uv})} be the number of edges.
#'   If selfLoops == F, then the code uses the approximation
#'   \deqn{
#'     Poisson(\lambda_{ij})
#'     \approx Binomial(M, \lambda_{ij}/\sum_{uv}\lambda_{uv})
#'   }
#'
#'   This approximation is good when total edges is ~n or larger and max
#'   \eqn{\lambda_{ij}} ~constant or smaller. this thresholds each edge;
#'     multiple edges are replaced by single edges.
#'
#' @references Rohe, Karl, Jun Tao, Xintian Han, and Norbert Binkiewicz. 2017.
#'    “A Note on Quickly Sampling a Sparse Matrix with Low Rank Expectation.”
#'    Journal of Machine Learning Research; 19(77):1−13, 2018.
#'    <http://www.jmlr.org/papers/v19/17-128.html>
#'
#' @export
#'
#' @examples
#'
#' # TODO
#'
fastRG <- function(X, S, Y = NULL, avgDeg = NULL, simple = FALSE,
                   PoissonEdges = TRUE, directed = TRUE,
                   selfLoops = TRUE, returnEdgeList = FALSE) {

  # input validation: X, S, Y must be non-negative
  if (any(X < 0) || any(S < 0) || any(Y < 0))
    stop("`X`, `S`, `Y` can only contain non-negative elements.", call. = FALSE)

  if (!directed & !is.null(Y))
    stop("Must not specify `Y` for undirected graphs.", call. = FALSE)

  # output will be asymmetric and potentially rectangular
  if (is.null(Y)) {
    Y <- X
  } else {
    directed <- TRUE
    selfLoops <- TRUE
    simple <- FALSE
  }

  n <- nrow(X)
  d <- nrow(Y)

  # more input validation

  # TODO: issue: d may not be defined here if Y is NULL, as is the default
  if (!directed && n != d)
    stop(
      "`n = nrow(X)` and `d = nrow(Y)` must be the same for undirected graphs",
      call. = FALSE
    )

  if (simple) {
    selfLoops <- FALSE
    directed <- FALSE
    PoissonEdges <- FALSE
  }

  # if avgDeg is specified, then scale S by the appropriate amount.
  if (!is.null(avgDeg)) {

    # find the expected average degree in the poisson graph
    # TODO: should S be symmetrized before this?
    # TODO: avgDeg vs avDeg is some confusing naming
    eDbar <- howManyEdges(X, S, Y)[["avDeg"]]

    S <- S * avgDeg / eDbar
  }

  # if undirected, symmetrize by setting S := (S + t(S))/2
  # then divide result by 2 because this doubles edge probabilities
  if (!directed)
    S <- (S + t(S)) / 4

  K1 <- ncol(X)
  K2 <- ncol(Y)

  Cx <- diag(colSums(X), nrow = K1, ncol = K1)
  Cy <- diag(colSums(Y), nrow = K2, ncol = K2)

  St <- Cx %*% S %*% Cy

  # number of edges to sample
  m <- stats::rpois(n = 1, lambda = sum(St))

  # if no edges, return empty matrix.
  if (m == 0) {

    if (returnEdgeList)
      return(matrix(NA, nrow = 0, ncol = 2))

    return(sparseMatrix(1:n, 1:d, x = 0, dims = c(n, d)))
  }

  # this simulates \varpi, denoted here as tabUV.
  # element u,v is the number of edges between column u and column v.

  tabUV <- matrix(stats::rmultinom(n = 1, size = m, prob = St), nrow = K1, ncol = K2)
  cumsumUV <- matrix(cumsum(tabUV), nrow = K1, ncol = K2)

  # cbind(eo, ei) is going to be the edge list
  #   eo: "edge out node".
  #   ei: "edge in node"

  eo <- rep(NA, m)
  ei <- eo

  # to avoid doing K1*K2 samples for I and J, we can instead take
  #   only K1 + K2 samples.  This requires some awkward indexing.
  #   for(u in 1:K1) for(v in 1:K2) tabUV[u,v]  <-
  #      eventually this is going to happen.
  #     because this first sets u = 1 and loops through
  #     the different values of v,
  #     the awkward indexing will only apply to the "v" or the ei vector.
  #     for eo, we can just loop through u in 1:K1...

  ticker <- 1
  blockDegreesU <- rowSums(tabUV)
  for (u in 1:K1) {
    if (blockDegreesU[u] > 0) {
      eo[ticker:(ticker + blockDegreesU[u] - 1)] <-
        sample(n, size = blockDegreesU[u], replace = T, prob = X[, u])
      ticker <- ticker + blockDegreesU[u]
    }
  }

  # for ei, things are more awkward.  if we had pointers, perhaps there
  # would be a faster way... instead create ei-tmp.  which will hold the
  # values of ei, but in the wrong order.  another loop will take care of
  # the indexing.

  eitmp <- ei
  ticker <- 1
  blockDegreesV <- colSums(tabUV)

  for (v in 1:K2) {
    if (blockDegreesV[v] > 0) {
      eitmp[ticker:(ticker + blockDegreesV[v] - 1)] <- sample(d, size = blockDegreesV[v], replace = T, prob = Y[, v])
      ticker <- ticker + blockDegreesV[v]
    }
  }

  # put ei-tmp in the correct order, to match up with eo
  ticker <- 1
  tickerV <- c(1, cumsum(blockDegreesV))
  for (u in 1:K1) {
    for (v in 1:K2) {
      edgesInThisBlock <- tabUV[u, v]
      if (edgesInThisBlock > 0) {
        ei[ticker:(ticker + edgesInThisBlock - 1)] <- eitmp[tickerV[v]:(tickerV[v] + edgesInThisBlock - 1)]
        tickerV[v] <- tickerV[v] + edgesInThisBlock
        ticker <- ticker + edgesInThisBlock
      }
    }
  }

  if (!selfLoops) {
    edges_to_self <- eo == ei
    eo <- eo[-edges_to_self]
    ei <- ei[-edges_to_self]
  }

  # symmetrize the edge list. this doubles the edge probabilities!
  if (!directed) {
    eo_copy <- eo
    eo <- c(eo, ei)
    ei <- c(ei, eo_copy)
  }

  if (returnEdgeList)
    return(cbind(eo, ei))

  if (PoissonEdges)
    # NOTE: x = 1 is correct to create a multigraph adjacency matrix
    # here. see ?Matrix::sparseMatrix for details, in particular the
    # documentation for arguments `i, j`
    A <- sparseMatrix(eo, ei, x = 1, dims = c(n, d), symmetric = !directed)
  else
    A <- sparseMatrix(eo, ei, dims = c(n, d), symmetric = !directed)

  A
}
