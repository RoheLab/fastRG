#' Sample a Random Dot Product Graph (RDPG)
#'
#' @param X An `n` by `k_1` matrix.
#' @param S A `k_1` by `k_2` matrix.
#' @param Y A `d` by `k_2` matrix. Defaults to `X`.
#'
#' @param avg_deg When specified, rescales parameter such that the
#'   expected degree is `avg_deg` in the Poisson multi-graph. When
#'   `poisson_edges = FALSE`, the resulting graph will have lower a
#'   average degree than `avg_deg` due to lack of multi-edges. When
#'   the graph is sparse, the expected number of edges for the Poisson
#'   multi-graph and Bernoulli graph are nearly the
#'   same. Defaults to `NULL`, such that no scaling occurs.
#'
#' @param simple When `TRUE`, indicates that you want to work with undirected
#'   graphs where self-loops and multi-edges are prohibited. Accomplishes
#'   this by setting `directed = FALSE`, `allow_self_loops = FALSE`, and
#'   `poisson_edges = FALSE`, and then ignoring arguments `directed`,
#'   `allow_self_loops` and `poisson_edges`. Defaults to `FALSE`.
#'
#' @param poisson_edges Logical indicating whether or not multi-edges are
#'   allowed. Defaults to `TRUE`, which keeps multi-edges and produces
#'   a multi-graph. When `FALSE`, only single edges are allowed, resulting
#'   in a graph. See **details** for some additional comments. Effected by
#'   `simple` argument.
#'
#' @param directed Logical indicating whether or not the graph should be
#'   directed. Default is `directed = TRUE`. When `directed = FALSE`,
#'   symmetrizes `S` internally. When `X = Y` (which is the default when
#'   no `Y` is specified), this results in a symmetric adjacency matrix
#'   as output. When `avg_deg` is specified and the desired graph is directed,
#'   the average degree scaling is on the out-degree of each node (or the
#'   row sums if you prefer to think in terms of the adjacency matrix).
#'   Effected by the `simple` argument.
#'
#' @param allow_self_loops Logical indicating whether edges are allowed from
#'   a node back to itself. Defaults to `TRUE`. When `FALSE`, sampling
#'   proceeds normally, and then self-loops are removed after sampling
#'   is completed. Effected by the `simple` argument.
#'
#' @param return_edge_list Logical indicating whether to return an edgelist
#'   rather than an adjacency matrix. Defaults to `FALSE`.
#'
#' @return A random Poisson (or Bernoulli) dot product graph. By default,
#'   returns a [Matrix::sparseMatrix()] in CSC format (i.e. of abstract
#'   class `CsparseMatrix`). When the graph is undirected, the `sparseMatrix`
#'   will also be symmetric.
#'
#'   If `return_edge_list = TRUE`, instead returns the graph as an edgelist
#'   in matrix format, with contains two columns called `from` and `to`.
#'   The `from` column contains the index of the outgoing edge,
#'   and the `to` column contains the index of the incoming edge.
#'
#' @details TODO: clean up the following:
#'
#'   This uses a poisson approximation to the binomial...
#'
#'   Let \eqn{M ~ Poisson(\sum_{uv} \lambda_{uv})} be the number of edges.
#'   If allow_self_loops == F, then the code uses the approximation
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
#'    Journal of Machine Learning Research; 19(77):1-13, 2018.
#'    <http://www.jmlr.org/papers/v19/17-128.html>
#'
#' @export
#'
#' @examples
#'
#' set.seed(372)
#'
#' n <- 10000
#' d <- 1000
#'
#' k1 <- 5
#' k2 <- 3
#'
#' X <- matrix(rpois(n = n * k1, 1), nrow = n)
#' Y <- matrix(rpois(n = d * k2, 1), nrow = d)
#' S <- matrix(runif(n = k1 * k2, 0, .1), nrow = K1)
#'
#' # all the bells and whistles, return graph as a matrix
#' A <- fastRG(X, S, Y, avg_deg = 10)
#'
#' # a nice binary, symmetric graph with no self-loops
#' B <- fastRG(X, S, avg_deg = 10, simple = TRUE)
#'
#' # get the graph as an edge list
#' edge_list <- fastRG(X, S, return_edge_list = TRUE)
#'
#' # a symmetric graph
#' C <- fastRG(X, S, directed = FALSE)
#'
fastRG <- function(X, S, Y = NULL, avg_deg = NULL, simple = FALSE,
                   poisson_edges = TRUE, directed = TRUE,
                   allow_self_loops = TRUE, return_edge_list = FALSE) {

  # TODO: check dimensions. also check in expected()

  if (any(X < 0) || any(S < 0) || any(Y < 0)) {
    stop("`X`, `S`, `Y` can only contain non-negative elements.", call. = FALSE)
  }

  if (!directed & !is.null(Y)) {
    stop("Must not specify `Y` for undirected graphs.", call. = FALSE)
  }

  # output will be asymmetric and potentially rectangular
  if (is.null(Y)) {
    Y <- X
  } else {
    directed <- TRUE
    allow_self_loops <- TRUE
    simple <- FALSE
  }

  n <- nrow(X)
  d <- nrow(Y)

  # more input validation

  # TODO: issue: d may not be defined here if Y is NULL, as is the default
  if (!directed && n != d) {
    stop(
      "`n = nrow(X)` and `d = nrow(Y)` must be the same for undirected graphs",
      call. = FALSE
    )
  }

  if (simple) {
    allow_self_loops <- FALSE
    directed <- FALSE
    poisson_edges <- FALSE
  }


  # if avg_deg is specified, then scale S by the appropriate amount.
  if (!is.null(avg_deg)) {

    # find the expected average degree in the poisson graph
    # TODO: should S be symmetrized before this?
    # TODO: avg_deg vs avDeg is some confusing naming
    eDbar <- expected(X, S, Y)$degree

    S <- S * avg_deg / eDbar
  }

  # if undirected, symmetrize by setting S := (S + t(S))/2
  # then divide result by 2 because this doubles edge probabilities
  if (!directed) {
    print(class(S))
    S <- (S + t(S)) / 4
  }

  K1 <- ncol(X)
  K2 <- ncol(Y)

  Cx <- diag(colSums(X), nrow = K1, ncol = K1)
  Cy <- diag(colSums(Y), nrow = K2, ncol = K2)

  St <- Cx %*% S %*% Cy

  # number of edges to sample
  m <- stats::rpois(n = 1, lambda = sum(St))

  # if no edges, return empty matrix.
  if (m == 0) {
    if (return_edge_list) {
      return(matrix(NA, nrow = 0, ncol = 2))
    }

    return(sparseMatrix(1:n, 1:d, x = 0, dims = c(n, d)))
  }

  # this simulates \varpi, denoted here as tabUV.
  # element u,v is the number of edges between column u and column v.

  tabUV <- matrix(stats::rmultinom(n = 1, size = m, prob = St), nrow = K1, ncol = K2)
  cumsumUV <- matrix(cumsum(tabUV), nrow = K1, ncol = K2)

  # cbind(eo, ei) is going to be the edge list
  #   eo: "edge out node".
  #   ei: "edge in node"

  eo <- integer(m)
  ei <- integer(m)

  print("boom")

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

  if (!allow_self_loops) {
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

  if (return_edge_list) {
    el <- cbind(eo, ei)
    colnames(el) <- c("from", "to")
    return(el)
  }

  if (poisson_edges) {
    # NOTE: x = 1 is correct to create a multigraph adjacency matrix
    # here. see ?Matrix::sparseMatrix for details, in particular the
    # documentation for arguments `i, j`
    A <- sparseMatrix(eo, ei, x = 1, dims = c(n, d), symmetric = !directed)
  } else {
    A <- sparseMatrix(eo, ei, dims = c(n, d), symmetric = !directed)
  }

  A
}
