new_undirected_overlapping_sbm <- function(
    X, S,
    Z,
    pi,
    B,
    sorted,
    ...,
    subclass = character()) {
  subclass <- c(subclass, "undirected_overlapping_sbm")
  overlapping_sbm <- undirected_factor_model(X, S, ..., subclass = subclass)
  overlapping_sbm$Z <- Z
  overlapping_sbm$pi <- pi
  overlapping_sbm$B <- B
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

  if (!is_invertible(values$S)) {
    stop(
      "`B` must be an invertible matrix.",
      call. = FALSE
    )
  }

  if (max(values$B) > 1 || min(values$B) < 0) {
    stop(
      "All elements of `B` must be contained in [0, 1].",
      call. = FALSE
    )
  }

  x
}

#' Create an undirected overlapping degree corrected stochastic blockmodel object
#'
#' To specify a overlapping stochastic blockmodel, you must specify
#' the number of nodes (via `n`),
#' the mixing matrix (via `k` or `B`), and the  block
#' probabilities (optional, via `pi`). We provide defaults for most of these
#' options to enable rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n The number of nodes in the overlapping SBM.
#'
#' @param k (mixing matrix) The number of blocks in the blockmodel.
#'   Use when you don't want to specify the mixing-matrix by hand.
#'   When `k` is specified, `B` is set to a diagonal dominant matrix with
#'   value `0.8` along the diagonal and `0.1 / (k - 1)` on the
#'   off-diagonal. `k` defaults to `NULL`. You must specify either `k`
#'   or `B`, but not both.
#'
#' @param B (mixing matrix) A `k` by `k` matrix of block connection
#'   probabilities. The probability that a node in block `i` connects
#'   to a node in community `j` is `Poisson(B[i, j])`. Must be
#'   an *invertible*, symmetric square matrix.
#'   `matrix` and `Matrix` objects are both
#'   acceptable. If `B` is not symmetric, it will be
#'   symmetrized via the update `B := B + t(B)`. Defaults to `NULL`.
#'   You must specify either `k` or `B`, but not both.
#'
#' @param pi (block probabilities) Probability of membership in each
#'   block. Membership in each block is independent under the
#'   overlapping SBM. Defaults to `rep(1 / k, k)`.
#'
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block. Useful for plotting.
#'   Defaults to `TRUE`.
#'
#' @param force_pure Logical indicating whether or not to force presence of
#'   "pure nodes" (nodes that belong only to a single community) for the sake
#'   of identifiability. To include pure nodes, block membership sampling
#'   first proceeds as per usual. Then, after it is complete, `k` nodes
#'   are chosen randomly as pure nodes, one for each block. Defaults to `TRUE`.
#'
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#'
#' @return An `undirected_overlapping_sbm` S3 object, a subclass of the
#'   [undirected_factor_model()] with the following additional
#'   fields:
#'
#'   - `pi`: Sampling probabilities for each block.
#'
#'   - `sorted`: Logical indicating where nodes are arranged by
#'     block (and additionally by degree heterogeneity parameter)
#'     within each block.
#'
#' @references
#'
#' Kaufmann, Emilie, Thomas Bonald, and Marc Lelarge.
#' "A Spectral Algorithm with Additive Clustering for the Recovery of
#' Overlapping Communities in Networks," Vol. 9925.
#' Lecture Notes in Computer Science.
#' Cham: Springer International Publishing, 2016.
#' https://doi.org/10.1007/978-3-319-46379-7.
#'
#' Latouche, Pierre, Etienne Birmelé, and Christophe Ambroise.
#' "Overlapping Stochastic Block Models with Application to the
#' French Political Blogosphere." The Annals of Applied Statistics 5,
#' no. 1 (March 2011): 309–36. https://doi.org/10.1214/10-AOAS382.
#'
#' Zhang, Yuan, Elizaveta Levina, and Ji Zhu. "Detecting
#' Overlapping Communities in Networks Using Spectral Methods."
#'  ArXiv:1412.3432, December 10, 2014. http://arxiv.org/abs/1412.3432.
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
#' ## Identifiability
#'
#' In order to be identifiable, an overlapping SBM must satisfy two conditions:
#'
#' 1. `B` must be invertible, and
#' 2. the must be at least one "pure node" in each block that belongs to no
#'   other blocks.
#'
#' ## Block memberships
#'
#' Note that some nodes may not belong to any blocks.
#'
#' **TODO**
#'
#' ## Edge formulation
#'
#' Once we know the block memberships, we need one more
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
#' k <- 5
#' n <- 1000
#' B <- matrix(stats::runif(k * k), nrow = k, ncol = k)
#'
#' pi <- c(1, 2, 4, 1, 1) / 5
#'
#' custom_overlapping_sbm <- overlapping_sbm(
#'   n = 200,
#'   B = B,
#'   pi = pi,
#'   expected_degree = 5
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
    n, k = NULL, B = NULL,
    ...,
    pi = rep(1 / k, k),
    sort_nodes = TRUE,
    force_pure = TRUE) {

  ### mixing matrix

  if (is.null(k) && is.null(B)) {
    stop("Must specify either `k` or `B`.", call. = FALSE)
  } else if (is.null(B)) {

    if (k < 1) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    message(
      "Setting `B` to a matrix with value 0.8 on the diagonal and ",
      "0.1 / (k - 1) on the off-diagonal. This parameterization may change ",
      "in the future. Explicitly set `B` for reproducible results."
    )

    B <- matrix(0.1 / (k - 1), nrow = k, ncol = k)
    diag(B) <- 0.8

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

  if (k > 1) {
    B <- B[order(pi), ]
    B <- B[, order(pi)]
  }

  pi <- sort(pi) # do not normalize pi!

  # sample block memberships

  X <- Matrix(0, nrow = n, ncol = k)
  colnames(X) <- paste0("block", 1:k)

  for (i in 1:k) {
    ind <- stats::rbinom(n, 1, pi[i])
    X[, i] <- ind
  }

  if (force_pure) {
    pure_indices <- sample(n, k)

    for (i in 1:k) {
      X[pure_indices[i], ] <- 0
      X[pure_indices[i], i] <- 1
    }
  }

  if (sort_nodes) {
    X <- sort_by_all_columns(X)
  }

  overlapping_sbm <- new_undirected_overlapping_sbm(
    X = X,
    S = B,  # accepts B but transforms by symmetrizing and scaling internally
    B = B,
    Z = X,
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
  cat("Block memberships (Z):", dim_and_class(x$Z), "\n")
  cat("Block probabilities (pi):", dim_and_class(x$pi), "\n\n")
  cat("Block connection propensities (B):", dim_and_class(x$B), "\n\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
