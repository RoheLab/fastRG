validate_undirected_sbm <- function(x) {
  values <- unclass(x)

  if (!inherits(x, "undirected_sbm")) {
    stop(
      "Undirected SBMs must inherit \"undirected_sbm\" class!",
      call. = FALSE
    )
  }

  distinct_theta_by_block <- sapply(
    split(values$theta, values$z),
    function(x) length(unique(x))
  )

  if (any(distinct_theta_by_block) > 1) {
    stop(
      "`theta` must be equal for all nodes in a given block.",
      call. = FALSE
    )
  }

  x
}

#' Create an undirected stochastic blockmodel object
#'
#' To specify a stochastic blockmodel, you must specify
#' the number of nodes (via `n`), the mixing matrix (via `k` or `B`),
#' and the relative block probabilites (optional, via `pi`).
#' We provide defaults for most of these options to enable
#' rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n The number of nodes in the network. Must be
#'   a positive integer. This argument is required.
#'
#' @inherit dcsbm params details
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#'
#' @return An `undirected_sbm` S3 object, which is a subclass of the
#'   [dcsbm()] object.
#'
#' @export
#' @family stochastic block models
#' @family undirected graphs
#'
#' @details
#'
#' A stochastic block is equivalent to a degree-corrected
#' stochastic blockmodel where the degree heterogeneity parameters
#' have all been set equal to 1.
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_sbm <- sbm(n = 100, k = 5, expected_density = 0.01)
#' lazy_sbm
#'
#' # by default we get a multigraph (i.e. multiple edges are
#' # allowed between the same two nodes). using bernoulli edges
#' # will with an adjacency matrix with only zeroes and ones
#'
#' bernoulli_sbm <- sbm(
#'   n = 500,
#'   k = 30,
#'   poisson_edges = FALSE,
#'   expected_degree = 8
#' )
#'
#' bernoulli_sbm
#'
#' edgelist <- sample_edgelist(bernoulli_sbm)
#' edgelist
#'
#' A <- sample_sparse(bernoulli_sbm)
#'
#' # only zeroes and ones!
#' sign(A)
#'
#'
#' # sbm with repeated eigenvalues
#'
#' # block sizes equal by default, needed to prevent variation in spectrum
#' # from variation in block sizes. also need B to have a single repeated
#' # eigenvalue
#'
#' repeated_eigen <- sbm(
#'   n = 100,
#'   B = diag(rep(0.8, 5)),
#'   expected_degree = 10
#' )
#'
#' # exactly repeated eigenvalues in the population
#' e <- eigs_sym(repeated_eigen)
#' e$values
#'
sbm <- function(
    n,
    k = NULL, B = NULL,
    ...,
    block_sizes = NULL,
    pi = NULL,
    sort_nodes = TRUE,
    poisson_edges = TRUE,
    allow_self_loops = TRUE) {
  sbm <- dcsbm(
    theta = rep(1, n),
    k = k,
    B = B,
    ...,
    block_sizes = block_sizes,
    pi = pi,
    sort_nodes = sort_nodes,
    force_identifiability = FALSE,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops,
    subclass = "undirected_sbm"
  )

  validate_undirected_sbm(sbm)
}

#' @method print undirected_sbm
#' @export
print.undirected_sbm <- function(x, ...) {
  cat(glue("Undirected Stochastic Blockmodel\n", .trim = FALSE))
  cat(glue("--------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k}\n\n", .trim = FALSE))

  cat("Traditional SBM parameterization:\n\n")
  cat("Block memberships (z):", dim_and_class(x$z), "\n")
  cat("Block probabilities (pi):", dim_and_class(x$pi), "\n")
  cat(glue("Edge distribution: {x$edge_distribution}\n\n", .trim = FALSE))

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat("Poisson edges:", as.character(x$poisson_edges), "\n")
  cat("Allow self loops:", as.character(x$allow_self_loops), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
