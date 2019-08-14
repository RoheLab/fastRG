#' Sample a Chung-Lu graph
#'
#' @param theta A vector containing expected degree of each node. That is,
#'   the resulting graph will have `length(theta)` nodes. Note that as
#'   `sum(theta)` approaches `length(theta)`, you will start to sample
#'   very dense graphs.
#'
#' @inheritDotParams fastRG
#' @inherit fastRG return
#'
#' @return Never returns Poisson edges.
#'
#' @export
#'
#' @seealso [fastRG()]
#' @family bernoulli graphs
#'
#' @examples
#'
#' set.seed(27)
#'
#' expected_degree <- c(0, 1, 0, 2, 0)
#'
#' A <- chung_lu(expected_degree)
#'
#' # out degree
#' rowSums(A)
#'
#' # get the random dot product model parameters
#' params <- chung_lu_params(expected_degree)
#'
chung_lu <- function(theta, ...) {
  params <- chung_lu_params(theta)
  fastRG(params$X, params$S, poisson_edges = FALSE, ...)
}

#' @rdname chung_lu
#' @export
chung_lu_params <- function(theta) {
  if (!is.vector(theta)) {
    stop("`theta` must be a vector.", call. = FALSE)
  }

  if (any(theta < 0)) {
    stop("Elements of `theta` must be non-negative.", call. = FALSE)
  }

  n <- length(theta)
  X <- matrix(theta, nrow = n, ncol = 1)
  S <- matrix(1, nrow = 1, ncol = 1)

  list(X = X, S = S, Y = X)
}
