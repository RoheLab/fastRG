new_directed_erdos_renyi <- function(
  X,
  S,
  Y,
  p,
  poisson_edges,
  allow_self_loops,
  ...
) {
  er <- directed_factor_model(
    X,
    S,
    Y,
    ...,
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

#' Create an directed erdos renyi object
#'
#' @param n Number of nodes in graph.
#'
#' @param p Probability of an edge between any two nodes, strictly
#'   between zero and one. You must specify either `p`,
#'   `expected_in_degree`, `expected_out_degree` or `expected_density`.
#'   When one of those is given, the stored `p` is the edge probability
#'   (the edge rate, for Poisson graphs) of the rescaled model, not the `p`
#'   given here.
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
  n,
  ...,
  p = NULL,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
) {
  X <- Matrix(1, nrow = n, ncol = 1)
  Y <- Matrix(1, nrow = n, ncol = 1)

  # the scaling arguments arrive via `...`, so look for them by name (a
  # bare `expected_in_degree` would find the generic function)
  scaling <- c("expected_in_degree", "expected_out_degree", "expected_density")

  if (is.null(p) && !any(scaling %in% names(list(...)))) {
    stop(
      "Must specify either `p`, `expected_in_degree`, ",
      "`expected_out_degree` or `expected_density`.",
      call. = FALSE
    )
  }

  if (is.null(p)) {
    p <- 0.5 # doesn't matter, will get rescaled anyway
  } else if (!is.numeric(p) || length(p) != 1 || p <= 0 || 1 <= p) {
    stop("`p` must be strictly between zero and one.", call. = FALSE)
  }

  S <- matrix(p, nrow = 1, ncol = 1)

  er <- new_directed_erdos_renyi(
    X,
    S,
    Y,
    p = p,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops,
    ...
  )

  # store the edge probability of the final, rescaled model
  if (is_rescaled(...)) {
    er$p <- edge_scale_S(er)[1, 1]
  }

  validate_directed_erdos_renyi(er)
}
