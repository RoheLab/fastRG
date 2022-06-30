new_undirected_erdos_renyi <- function(X, S, p, ...) {
  er <- undirected_factor_model(
    X, S, ...,
    subclass = "undirected_erdos_renyi"
  )
  er$p <- p
  er
}

validate_undirected_erdos_renyi <- function(x) {

  values <- unclass(x)

  if (ncol(values$X) != 1) {
    stop("`X` must have a single column.", call. = FALSE)
  }

  if (values$p <= 0 || 1 <= values$p) {
    stop("`p` must be strictly between zero and one.", call. = FALSE)
  }

  x
}

#' Create an undirected erdos renyi object
#'
#' @param n Number of nodes in graph.
#'
#' @param p Probability of an edge between any two nodes. You must specify
#'   either `p` or `expected_degree`.
#'
#' @inheritDotParams undirected_factor_model expected_degree
#' @inherit undirected_factor_model params return
#'
#' @export
#'
#' @family erdos renyi
#' @family undirected graphs
#'
#' @examples
#'
#' set.seed(87)
#'
#' er <- erdos_renyi(n = 10, p = 0.1)
#' er
#'
#'
#' er <- erdos_renyi(n = 10, expected_density = 0.1)
#' er
#'
#' big_er <- erdos_renyi(n = 10^6, expected_degree = 5)
#' big_er
#'
#' A <- sample_sparse(er)
#' A
#'
erdos_renyi <- function(
  n, ..., p = NULL,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  X <- Matrix(1, nrow = n, ncol = 1)

  if (is.null(p) && is.null(expected_degree)) {
    stop("Must specify either `p` or `expected_degree`.", call. = FALSE)
  }

  if (is.null(p)) {
    p <- 0.5  # doesn't matter, will get rescaled anyway
  }

  B <- matrix(p)

  er <- new_undirected_erdos_renyi(X, B, p = p, ...,
                                   poisson_edges = poisson_edges,
                                   allow_self_loops = allow_self_loops)
  validate_undirected_erdos_renyi(er)
}
