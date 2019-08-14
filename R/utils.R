#' Calculate the expected number of edges in Poisson RDPG graph
#'
#' @param X An `n` by `k_1` matrix.
#' @param S A `k_1` by `k_2` matrix.
#' @param Y A `d` by `k_2` matrix. Defaults to `X`.
#'
#' @return A named list with three measurements from a Poisson RDPG(X, S, Y)
#'     multi-graph:
#'
#'   - `count`: Expected number of edges.
#'   - `degree`: Expected degree of each node.
#'   - `density`: Expected edge density. Equals the expected number of edges
#'     divided by the total number of edges.
#'
#' @details Note that the running time of [fastRG()] is proportional to
#'   `count`, the expected number of edges in the graph.
#'
#' @export
#'
#' @examples
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
#' expected(X, S, Y)
expected <- function(X, S, Y = X) {

  print("start of expected")

  if (any(X < 0) || any(S < 0) || any(Y < 0)) {
    stop("`X`, `S`, `Y` can only contain non-negative elements.", call. = FALSE)
  }

  Cx <- diag(colSums(X), nrow = ncol(X), ncol = ncol(X))
  Cy <- diag(colSums(Y), nrow = ncol(Y), ncol = ncol(Y))

  # convert integers to numerics to fend off integer overflow
  nx <- as.numeric(nrow(X))
  ny <- as.numeric(nrow(Y))

  count <- sum(Cx %*% S %*% Cy)
  degree <- count / nx

  # TODO: set `density = NULL` and issue a warning on overflows

  print("end of expected")

  list(
    count = count,
    degree = degree,
    density = count / (nx * ny)
  )
}
