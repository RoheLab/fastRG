#' Create an undirected Chung-Lu object
#'
#' To specify a Chung-Lu graph, you must specify
#' the degree-heterogeneity parameters (via `n` or `theta`).
#' We provide reasonable defaults to enable rapid exploration
#' or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n (degree heterogeneity) The number of nodes in the graph.
#'   Use when you don't want to specify the degree-heterogeneity
#'   parameters `theta` by hand. When `n` is specified, `theta`
#'   is randomly generated from a `LogNormal(2, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `n` defaults to `NULL`. You must specify either `n`
#'   or `theta`, but not both.
#'
#' @param theta (degree heterogeneity) A numeric vector
#'   explicitly specifying the degree heterogeneity
#'   parameters. This implicitly determines the number of nodes
#'   in the resulting graph, i.e. it will have `length(theta)` nodes.
#'   Must be positive. Setting to a vector of ones recovers
#'   an erdos renyi graph.
#'   Defaults to `NULL`. You must specify either `n` or `theta`,
#'   but not both.
#'
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block and by `theta`. Useful for plotting.
#'   Defaults to `TRUE`.
#'
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#' @inheritParams undirected_factor_model
#' @inheritParams dcsbm
#'
#' @return An `undirected_chung_lu` S3 object, a subclass of [dcsbm()].
#'
#' @export
#' @family undirected graphs
#'
#' @examples
#'
#' set.seed(27)
#'
#' cl <- chung_lu(n = 1000, expected_density = 0.01)
#' cl
#'
#' theta <- round(stats::rlnorm(100, 2))
#'
#' cl2 <- chung_lu(
#'   theta = theta,
#'   expected_degree = 5
#' )
#'
#' cl2
#'
#' edgelist <- sample_edgelist(cl)
#' edgelist
#'
chung_lu <- function(
  n = NULL, theta = NULL,
  ...,
  sort_nodes = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE,
  force_identifiability = FALSE) {

  ### degree heterogeneity parameters

  if (is.null(n) && is.null(theta)) {
    stop("Must specify either `n` or `theta`.", call. = FALSE)
  } else if (is.null(theta)) {

    if (n < 1) {
      stop("`n` must be a positive integer.", call. = FALSE)
    }

    message(
      "Generating random degree heterogeneity parameters `theta` from a ",
      "LogNormal(2, 1) distribution. This distribution may change ",
      "in the future. Explicitly set `theta` for reproducible results.\n"
    )

    theta <- stats::rlnorm(n, meanlog = 2, sdlog = 1)
  } else if (is.null(n)) {
    n <- length(theta)
  }

  chung_lu <- dcsbm(
    B = matrix(1),
    theta = theta,
    pi = 1,
    sort_nodes = sort_nodes,
    poisson_edges = poisson_edges,
    force_identifiability = force_identifiability,
    allow_self_loops = allow_self_loops,
    ...
  )

  # don't validate the chung-ln-ness because there is no low level
  # chung-lu constructor
  validate_undirected_dcsbm(chung_lu)
}

#' @method print undirected_chung_lu
#' @export
print.undirected_chung_lu <- function(x, ...) {

  cat(glue("Undirected Chung-Lu Graph\n", .trim = FALSE))
  cat(glue("-------------------------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))

  cat("Traditional Chung-Lu parameterization:\n\n")
  cat("Degree heterogeneity (theta):", dim_and_class(x$theta), "\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat("Poisson edges:", as.character(x$poisson_edges), "\n")
  cat("Allow self loops:", as.character(x$allow_self_loops), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
