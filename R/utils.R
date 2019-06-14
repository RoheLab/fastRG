#' Calculate the expected number of edges in Poisson gRDPG graph
#'
#' @param X An `n` by `k_1` matrix.
#' @param S A `k_1` by `k_2` matrix.
#' @param Y A `d` by `k_2` matrix. Defaults to `X`.
#'
#' @return A named list with three measurements from a Poisson gRDPG(X, S, Y)
#'     multi-graph:
#'   - `em`: Expected number of edges.
#'   - `avgDeg`: Expected degree of each node.
#'   - `density`: Expected edge density. Equals the expected number of edges
#'     divided by the total number of edges.
#'
#' @details Note that the running time of [fastRG()] is proportional to the
#'   expected number of edges `em`.
#' @export
#'
#' @examples
#'
#' # TODO
#'
howManyEdges <- function(X, S, Y = X) {

  Cx <- diag(colSums(X), nrow = ncol(X), ncol = ncol(X))
  Cy <- diag(colSums(Y), nrow = ncol(Y), ncol = ncol(Y))

  # convert integers to numerics to fend off integer overflow
  nx <- as.numeric(nrow(X))
  ny <- as.numeric(nrow(Y))

  em <- sum(Cx %*% S %*% Cy)
  avDeg <- em / nx

  # TODO: set `density = NULL` and issue a warning on overflows

  list(
    em = em,
    avDeg = avDeg,
    density = em / (nx * ny)
  )
}
