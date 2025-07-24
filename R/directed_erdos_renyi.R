new_directed_erdos_renyi <- function(X, S, Y, p, poisson_edges, allow_self_loops, ...) {
  er <- directed_factor_model(
    X, S, Y, ...,
    subclass = "directed_erdos_renyi",
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )

  er$p <- p
  er
}

validate_directed_erdos_renyi <- function(x) {
  values <- unclass(x)

  if (ncol(values$X) != 1) {
    stop("`X` must have a single column.", call. = FALSE)
  }

  if (ncol(values$Y) != 1) {
    stop("`Y` must have a single column.", call. = FALSE)
  }

  if (values$p < 0) {
    stop("`p` must be strictly non-negative.", call. = FALSE)
  }

  x
}

#' Create an directed erdos renyi object
#'
#' @param n Number of nodes in graph.
#'
#' @param p Probability of an edge between any two nodes. You must specify
#'   either `p`, `expected_in_degree`, or `expected_out_degree`.
#'
#' @inheritDotParams directed_factor_model expected_in_degree expected_out_degree
#' @inherit directed_factor_model params return
#'
#' @export
#'
#' @family erdos renyi
#' @family directed graphs
#'
#' @examples
#'
#' set.seed(87)
#'
#' er <- directed_erdos_renyi(n = 10, p = 0.1)
#' er
#'
#' big_er <- directed_erdos_renyi(n = 1000, expected_in_degree = 5)
#' big_er
#'
#' A <- sample_sparse(er)
#' A
#'
directed_erdos_renyi <- function(
    n, ..., p = NULL,
    poisson_edges = TRUE,
    allow_self_loops = TRUE) {
  X <- Matrix(1, nrow = n, ncol = 1)
  Y <- Matrix(1, nrow = n, ncol = 1)

  if (is.null(p) && is.null(expected_in_degree) &&
    is.null(expected_out_degree)) {
    stop(
      "Must specify either `p`, `expected_in_degree` or ",
      " `expected_out_degree`.",
      call. = FALSE
    )
  }

  if (is.null(p)) {
    p <- 0.5 # doesn't matter, will get rescaled anyway
  }

  S <- matrix(p, nrow = 1, ncol = 1)

  er <- new_directed_erdos_renyi(
    X, S, Y,
    p = p,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops,
    ...
  )

  validate_directed_erdos_renyi(er)
}
