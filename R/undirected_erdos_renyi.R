new_undirected_erdos_renyi <- function(X, S, p, ...) {
  er <- undirected_factor_model(
    X,
    S,
    ...,
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

  if (!is.numeric(values$p) || length(values$p) != 1 || values$p <= 0) {
    stop("`p` must be a positive scalar.", call. = FALSE)
  }

  if (!values$poisson_edges && values$p >= 1) {
    stop(
      "`p` must be strictly between zero and one for Bernoulli graphs.",
      call. = FALSE
    )
  }

  x
}

#' Create an undirected erdos renyi object
#'
#' @param n Number of nodes in graph.
#'
#' @param p Probability of an edge between any two nodes, strictly
#'   between zero and one. You must specify either `p`, `expected_degree`
#'   or `expected_density`. When `expected_degree` or `expected_density` is
#'   given, the stored `p` is the edge probability (the edge rate, for
#'   Poisson graphs) of the rescaled model, not the `p` given here.
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
#' big_er <- erdos_renyi(n = 1000, expected_degree = 5)
#' big_er
#'
#' A <- sample_sparse(er)
#' A
#'
erdos_renyi <- function(
  n,
  ...,
  p = NULL,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
) {
  X <- Matrix(1, nrow = n, ncol = 1)

  # `expected_degree` and `expected_density` arrive via `...`, so look for
  # them by name (a bare `expected_degree` would find the generic function)
  scaling <- c("expected_degree", "expected_density")

  if (is.null(p) && !any(scaling %in% names(list(...)))) {
    stop(
      "Must specify either `p`, `expected_degree` or `expected_density`.",
      call. = FALSE
    )
  }

  if (is.null(p)) {
    p <- 0.5 # doesn't matter, will get rescaled anyway
  } else if (!is.numeric(p) || length(p) != 1 || p <= 0 || 1 <= p) {
    stop("`p` must be strictly between zero and one.", call. = FALSE)
  }

  B <- matrix(p)

  er <- new_undirected_erdos_renyi(
    X,
    B,
    p = p,
    ...,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )

  # store the edge probability of the final, rescaled model
  if (is_rescaled(...)) {
    er$p <- edge_scale_S(er)[1, 1]
  }

  validate_undirected_erdos_renyi(er)
}
