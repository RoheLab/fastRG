new_undirected_overlapping_sbm <- function(
    X, S,
    theta,
    pi,
    sorted,
    ...,
    subclass = character()) {
  subclass <- c(subclass, "undirected_overlapping_sbm")
  overlapping_sbm <- undirected_factor_model(X, S, ..., subclass = subclass)
  overlapping_sbm$theta <- theta
  overlapping_sbm$pi <- pi
  overlapping_sbm$sorted <- sorted
  overlapping_sbm
}

validate_undirected_overlapping_sbm <- function(x) {

  values <- unclass(x)

  if (any(values$pi < 0) || any(values$pi > 1)) {
    stop(
      "All elements of `pi` must be contained in [0, 1].",
      call. = FALSE
    )
  }

  # TODO: this probably should be done more thoroughly

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

  x
}

#' Create an undirected degree corrected stochastic blockmodel object
#'
#' To specify a degree-corrected stochastic blockmodel, you must specify
#' the degree-heterogeneity parameters (via `n` or `theta`),
#' the mixing matrix (via `k` or `B`), and the relative block
#' probabilites (optional, via `pi`). We provide sane defaults for most of these
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
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#'
#' @return An `undirected_overlapping_sbm` S3 object, a subclass of the
#'   [undirected_factor_model()] with the following additional
#'   fields:
#'
#'   - `theta`: A numeric vector of degree-heterogeneity parameters.
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
#' overlapping stochastic blockmodel. First, for each node, we
#' independently determine if that node is a member of each block. This is
#' handled by `overlapping_sbm()`. Then, given these block memberships,
#' we randomly sample edges between nodes. This second
#' operation is handled by [sample_edgelist()],
#' [sample_sparse()], [sample_igraph()] and
#' [sample_tidygraph()], depending depending on your desired
#' graph representation.
#'
#' ## Block memberships
#'
#' **TODO**
#'
#' ## Degree heterogeneity
#'
#' In addition to block membership, the overlapping SBM also allows
#' nodes to have different propensities for edge formation.
#' We represent this propensity for node \eqn{i} by a positive
#' number \eqn{\theta_i}.
#'
#' **TODO** identifiability constraint?
#'
#' ## Edge formulation
#'
#' Once we know the block memberships and the degree
#' heterogeneity parameters \eqn{theta}, we need one more
#' ingredient, which is the baseline intensity of connections
#' between nodes in block `i` and block `j`. Then each edge
#' \eqn{A_{i,j}} is Poisson distributed with parameter
#'
#' **TODO**
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_overlapping_sbm <- overlapping_sbm(n = 1000, k = 5, expected_density = 0.01)
#' lazy_overlapping_sbm
#'
#' # sometimes you gotta let the world burn and
#' # sample a wildly dense graph
#'
#' dense_lazy_overlapping_sbm <- overlapping_sbm(n = 500, k = 3, expected_density = 0.8)
#' dense_lazy_overlapping_sbm
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
#' pi <- c(1, 2, 4, 1, 1) / 5
#'
#' custom_overlapping_sbm <- overlapping_sbm(
#'   theta = theta,
#'   B = B,
#'   pi = pi,
#'   expected_degree = 50
#' )
#'
#' custom_overlapping_sbm
#'
#' edgelist <- sample_edgelist(custom_overlapping_sbm)
#' edgelist
#'
#' # efficient eigendecompostion that leverages low-rank structure in
#' # E(A) so that you don't have to form E(A) to find eigenvectors,
#' # as E(A) is typically dense. computation is
#' # handled via RSpectra
#'
#' population_eigs <- eigs_sym(custom_overlapping_sbm)
#'
overlapping_sbm <- function(
    n = NULL, theta = NULL,
    k = NULL, B = NULL,
    ...,
    pi = rep(1 / k, k),
    sort_nodes = TRUE) {

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

    B <- Matrix(data = stats::runif(k * k), nrow = k, ncol = k)

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

  if (any(pi < 0) || any(pi > 1)) {
    stop("All elements of `pi` must be contained in [0, 1].", call. = FALSE)
  }

  # order mixing matrix by expected group size

  B <- B[order(pi), ]
  B <- B[, order(pi)]
  pi <- sort(pi)        # do not normalize pi!

  # sample block memberships


  X <- Matrix(0, nrow = n, ncol = k)
  colnames(X) <- paste0("block", 1:k)

  for (i in 1:k) {
    X[, i] <- stats::rbinom(n, 1, pi[i]) * theta
  }

  if (sort_nodes) {

    args <- as.list(as.data.frame(as.matrix(X)))
    args$decreasing <- TRUE

    X <- X[do.call(order, args),]
  }

  overlapping_sbm <- new_undirected_overlapping_sbm(
    X = X,
    S = B,
    theta = theta,
    pi = pi,
    sorted = sort_nodes,
    ...
  )

  validate_undirected_overlapping_sbm(overlapping_sbm)
}

#' @method print undirected_overlapping_sbm
#' @export
print.undirected_overlapping_sbm <- function(x, ...) {

  cat(glue("Undirected Degree-Corrected Overlapping Blockmodel\n", .trim = FALSE))
  cat(glue("-------------------------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "sorted" else "sorted"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k}\n\n", .trim = FALSE))

  cat("Traditional Overlapping SBM parameterization:\n\n")
  cat("Degree heterogeneity (theta):", dim_and_class(x$theta), "\n")
  cat("Block probabilities (pi):", dim_and_class(x$pi), "\n\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
