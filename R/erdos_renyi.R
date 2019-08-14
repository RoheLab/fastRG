#' Sample an Erdos-Renyi graph
#'
#' @param n Number of nodes in graph.
#'
#' @param p Probability of an edge between any two nodes. You must specify
#'   either `p` or `avg_deg`. If you do not specify `p`, uses `avg_deg / n`.
#'
#' @param avg_deg Desired average out degree. When specified, rescales
#'   sampling probabilities to achieve the desired average out degree.
#'   Defaults to `NULL`, such that there is no rescaling.
#'
#' @param directed Defaults to `FALSE` for Erdos-Renyi graphs. The default
#'   in the more general [fastRG()] is `TRUE`.
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
#' # get the random dot product parameters
#' erdos_renyi_params(n = 10, p = 0.1)
#'
#' # sample a small graph
#' A <- erdos_renyi(n = 10, p = 0.1)
#' A
#'
#' # out degree of each node
#' rowSums(A)
#'
#' # in degree of each node
#' colSums(A)
#'
#' # sample a much larger graph
#' B <- erdos_renyi(n = 10^6, avg_deg = 5)
#'
erdos_renyi <- function(n, p = NULL, avg_deg = NULL, directed = FALSE, ...) {
  params <- erdos_renyi_params(
    n = n, p = p, avg_deg = avg_deg,
    directed = directed, ...
  )
  fastRG(params$X, params$S, poisson_edges = FALSE, directed = directed, ...)
}

#' @rdname erdos_renyi
#' @export
erdos_renyi_params <- function(n, p = NULL, avg_deg = NULL, directed = FALSE, ...) {
  X <- matrix(1, nrow = n, ncol = 1)

  if (is.null(p) && is.null(avg_deg)) {
    stop("Must specify either `avg_deg` or `p`.", call. = FALSE)
  }

  if (is.null(p)) {
    p <- avg_deg / n
  }

  poisson_p <- -log(1 - p)
  S <- matrix(poisson_p, nrow = 1, ncol = 1)

  list(X = X, S = S, Y = X)
}
