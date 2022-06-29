new_undirected_dcsbm <- function(
  X, S,
  theta,
  z,
  pi,
  sorted,
  ...,
  subclass = character()) {
  subclass <- c(subclass, "undirected_dcsbm")
  dcsbm <- undirected_factor_model(X, S, ..., subclass = subclass)
  dcsbm$theta <- theta
  dcsbm$z <- z
  dcsbm$pi <- pi
  dcsbm$sorted <- sorted
  dcsbm
}

validate_undirected_dcsbm <- function(x) {

  values <- unclass(x)

  if (!is.factor(values$z)) {
    stop("`z` must be a factor.", call. = FALSE)
  }

  if (length(values$z) != nrow(values$X)) {
    stop("There must be one element of `z` for each row of `X`.")
  }

  if (!is.numeric(values$theta)) {
    stop("`theta` must be a numeric.", call. = FALSE)
  }

  if (length(values$theta) != nrow(values$X)) {
    stop(
      "There must be one element of `theta` for each row of `X`.",
      call. = FALSE
    )
  }

  if (any(values$theta < 0)) {
    stop("Elements of `theta` must be strictly positive.", call. = FALSE)
  }

  if (length(levels(values$z)) != values$k) {
    stop(
      "The number of levels of `z` must match the rank of the model.",
      call. = FALSE
    )
  }

  if (length(levels(values$z)) != values$k) {
    stop(
    "The number of levels of `z` must match the rank of the model.",
    call. = FALSE
    )
  }

  if (min(values$S) < 0) {
    stop(
      "All elements of `B` must be non-negative.",
      call. = FALSE
    )
  }

  x
}

#' Create an undirected degree corrected stochastic blockmodel object
#'
#' To specify a degree-corrected stochastic blockmodel, you must specify
#' the degree-heterogeneity parameters (via `n` or `theta`),
#' the mixing matrix (via `k` or `B`), and the relative block
#' probabilities (optional, via `pi`). We provide defaults for most of these
#' options to enable rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n (degree heterogeneity) The number of nodes in the blockmodel.
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
#'   a stochastic blockmodel without degree correction.
#'   Defaults to `NULL`. You must specify either `n` or `theta`,
#'   but not both.
#'
#' @param k (mixing matrix) The number of blocks in the blockmodel.
#'   Use when you don't want to specify the mixing-matrix by hand.
#'   When `k` is specified, the elements of `B` are drawn
#'   randomly from a `Uniform(0, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `k` defaults to `NULL`. You must specify either `k`
#'   or `B`, but not both.
#'
#' @param B (mixing matrix) A `k` by `k` matrix of block connection
#'   probabilities. The probability that a node in block `i` connects
#'   to a node in community `j` is `Poisson(B[i, j])`. Must be
#'   a square matrix. `matrix` and `Matrix` objects are both
#'   acceptable. If `B` is not symmetric, it will be
#'   symmetrized via the update `B := B + t(B)`. Defaults to `NULL`.
#'   You must specify either `k` or `B`, but not both.
#'
#' @param pi (relative block probabilities) Relative block
#'   probabilities. Must be positive, but do not need to sum
#'   to one, as they will be normalized internally.
#'   Must match the dimensions of `B` or `k`. Defaults to
#'   `rep(1 / k, k)`, or a balanced blocks.
#'
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block and by `theta`. Useful for plotting.
#'   Defaults to `TRUE`.
#'
#' @param force_identifiability Logical indicating whether or not to
#'   normalize `theta` such that it sums to one within each block. Defaults
#'   to `TRUE`.
#'
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#'
#' @return An `undirected_dcsbm` S3 object, a subclass of the
#'   [undirected_factor_model()] with the following additional
#'   fields:
#'
#'   - `theta`: A numeric vector of degree-heterogeneity parameters.
#'
#'   - `z`: The community memberships of each node, as a [factor()].
#'     The factor will have `k` levels, where `k` is the number of
#'     communities in the stochastic blockmodel. There will not
#'     always necessarily be observed nodes in each community.
#'
#'   - `pi`: Sampling probabilities for each block.
#'
#'   - `sorted`: Logical indicating where nodes are arranged by
#'     block (and additionally by degree heterogeneity parameter)
#'     within each block.
#'
#' @export
#' @family stochastic block models
#' @family undirected graphs
#'
#' @details
#'
#' # Generative Model
#'
#' There are two levels of randomness in a degree-corrected
#' stochastic blockmodel. First, we randomly chose a block
#' membership for each node in the blockmodel. This is
#' handled by `dcsbm()`. Then, given these block memberships,
#' we randomly sample edges between nodes. This second
#' operation is handled by [sample_edgelist()],
#' [sample_sparse()], [sample_igraph()] and
#' [sample_tidygraph()], depending depending on your desired
#' graph representation.
#'
#' ## Block memberships
#'
#' Let \eqn{z_i} represent the block membership of node \eqn{i}.
#' To generate \eqn{z_i} we sample from a categorical
#' distribution (note that this is a special case of a
#' multinomial) with parameter \eqn{\pi}, such that
#' \eqn{\pi_i} represents the probability of ending up in
#' the ith block. Block memberships for each node are independent.
#'
#' ## Degree heterogeneity
#'
#' In addition to block membership, the DCSBM also allows
#' nodes to have different propensities for edge formation.
#' We represent this propensity for node \eqn{i} by a positive
#' number \eqn{\theta_i}. Typically the \eqn{\theta_i} are
#' constrained to sum to one for identifiability purposes,
#' but this doesn't really matter during sampling (i.e.
#' without the sum constraint scaling \eqn{B} and \eqn{\theta}
#' has the same effect on edge probabilities, but whether
#' \eqn{B} or \eqn{\theta} is responsible for this change
#' is uncertain).
#'
#' ## Edge formulation
#'
#' Once we know the block memberships \eqn{z} and the degree
#' heterogeneity parameters \eqn{theta}, we need one more
#' ingredient, which is the baseline intensity of connections
#' between nodes in block `i` and block `j`. Then each edge
#' \eqn{A_{i,j}} is Poisson distributed with parameter
#'
#' \deqn{
#'   \lambda[i, j] = \theta_i \cdot B_{z_i, z_j} \cdot \theta_j.
#' }{
#'   \lambda_{i, j} = \theta[i] * B[z[i], z[j]] * \theta[j].
#' }
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_dcsbm <- dcsbm(n = 1000, k = 5, expected_density = 0.01)
#' lazy_dcsbm
#'
#' # sometimes you gotta let the world burn and
#' # sample a wildly dense graph
#'
#' dense_lazy_dcsbm <- dcsbm(n = 500, k = 3, expected_density = 0.8)
#' dense_lazy_dcsbm
#'
#' # explicitly setting the degree heterogeneity parameter,
#' # mixing matrix, and relative community sizes rather
#' # than using randomly generated defaults
#'
#' k <- 5
#' n <- 1000
#' B <- matrix(stats::runif(k * k), nrow = k, ncol = k)
#'
#' theta <- round(stats::rlnorm(n, 2))
#'
#' pi <- c(1, 2, 4, 1, 1)
#'
#' custom_dcsbm <- dcsbm(
#'   theta = theta,
#'   B = B,
#'   pi = pi,
#'   expected_degree = 50
#' )
#'
#' custom_dcsbm
#'
#' edgelist <- sample_edgelist(custom_dcsbm)
#' edgelist
#'
#' # efficient eigendecompostion that leverages low-rank structure in
#' # E(A) so that you don't have to form E(A) to find eigenvectors,
#' # as E(A) is typically dense. computation is
#' # handled via RSpectra
#'
#' population_eigs <- eigs_sym(custom_dcsbm)
#'
dcsbm <- function(
  n = NULL, theta = NULL,
  k = NULL, B = NULL,
  ...,
  pi = rep(1 / k, k),
  sort_nodes = TRUE,
  force_identifiability = TRUE) {

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

  ### mixing matrix

  if (is.null(k) && is.null(B)) {
    stop("Must specify either `k` or `B`.", call. = FALSE)
  } else if (is.null(B)) {

    if (k < 1) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    message(
      "Generating random mixing matrix `B` with independent ",
      "Uniform(0, 1) entries. This distribution may change ",
      "in the future. Explicitly set `B` for reproducible results."
    )

    B <- matrix(data = stats::runif(k * k), nrow = k, ncol = k)

  } else if (is.null(k)) {

    if (nrow(B) != ncol(B)) {
      stop("`B` must be a square matrix.", call. = FALSE)
    }

    k <- nrow(B)
  }

  ### block membership

  if (length(pi) != nrow(B) || length(pi) != ncol(B)) {
    stop("Length of `pi` must match dimensions of `B`.", call. = FALSE)
  }

  if (any(pi < 0)) {
    stop("All elements of `pi` must be >= 0.", call. = FALSE)
  }

  # order mixing matrix by expected group size

  if (k > 1) {
    B <- B[order(pi), ]
    B <- B[, order(pi)]
  }

  pi <- sort(pi / sum(pi))

  # sample block memberships

  z <- sample(k, n, replace = TRUE, prob = pi)
  z <- factor(z, levels = 1:k, labels = paste0("block", 1:k))

  if (sort_nodes) {
    z <- sort(z)
  }

  if (k > 1) {
    X <- sparse.model.matrix(~z + 0)
  } else {
    X <- Matrix(1, nrow = n, ncol = 1)
  }

  if (force_identifiability) {
    theta <- l1_normalize_within(theta, z)
  }

  X@x <- theta

  if (sort_nodes) {
    # note that X and z indexing must match
    X <- sort_by_all_columns(X)
  }

  if (k > 1) {
    theta <- X@x
  } else {
    theta <- sort(theta, decreasing = TRUE)
  }

  dcsbm <- new_undirected_dcsbm(
    X = X,
    S = B,
    theta = theta,
    z = z,
    pi = pi,
    sorted = sort_nodes,
    ...
  )

  validate_undirected_dcsbm(dcsbm)
}

#' @method print undirected_dcsbm
#' @export
print.undirected_dcsbm <- function(x, ...) {

  cat(glue("Undirected Degree-Corrected Stochastic Blockmodel\n", .trim = FALSE))
  cat(glue("-------------------------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k}\n\n", .trim = FALSE))

  cat("Traditional DCSBM parameterization:\n\n")
  cat("Block memberships (z):", dim_and_class(x$z), "\n")
  cat("Degree heterogeneity (theta):", dim_and_class(x$theta), "\n")
  cat("Block probabilities (pi):", dim_and_class(x$pi), "\n\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
