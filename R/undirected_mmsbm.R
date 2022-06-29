validate_undirected_mmsbm <- function(x) {

  values <- unclass(x)

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


  if (min(values$S) < 0) {
    stop(
      "All elements of `S` must be non-negative.",
      call. = FALSE
    )
  }

  if (!isTRUE(all.equal(rowSums(values$Z), rep(1, values$n)))) {
    stop(
      "Rows of `Z` must all sum to one.",
      call. = FALSE
    )
  }

  x
}

new_undirected_mmsbm <- function(
    X, S,
    theta,
    Z,
    alpha,
    sorted,
    ...,
    subclass = character()) {

  subclass <- c(subclass, "undirected_mmsbm")
  mmsbm <- undirected_factor_model(X, S, ..., subclass = subclass)
  mmsbm$theta <- theta
  mmsbm$alpha <- alpha
  mmsbm$Z <- Z
  mmsbm$sorted <- sorted
  mmsbm
}
#' Create an undirected degree-corrected mixed membership stochastic blockmodel object
#'
#' To specify a degree-corrected mixed membership stochastic blockmodel, you must specify
#' the degree-heterogeneity parameters (via `n` or `theta`),
#' the mixing matrix (via `k` or `B`), and the relative block
#' propensities (optional, via `alpha`). We provide defaults for most of these
#' options to enable rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param alpha (relative block propensities) Relative block
#'   propensities, which are parameters of a Dirichlet distribution.
#'   All elments of `alpha` must thus be positive.
#'   Must match the dimensions of `B` or `k`. Defaults to
#'   `rep(1, k)`, or balanced membership across blocks.
#'
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#' @inheritParams dcsbm
#' @inheritParams overlapping_sbm
#'
#' @return An `undirected_mmsbm` S3 object, a subclass of the
#'   [undirected_factor_model()] with the following additional
#'   fields:
#'
#'   - `theta`: A numeric vector of degree-heterogeneity parameters.
#'
#'   - `Z`: The community memberships of each node, a [matrix()] with
#'     `k` columns, whose row sums all equal one.
#'
#'   - `alpha`: Community membership proportion propensities.
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
#' stochastic blockmodel. First, we randomly choose how much
#' each node belongs to each block in the blockmodel. Each node
#' is one unit of block membership to distribute. This is
#' handled by `mmsbm()`. Then, given these block memberships,
#' we randomly sample edges between nodes. This second
#' operation is handled by [sample_edgelist()],
#' [sample_sparse()], [sample_igraph()] and
#' [sample_tidygraph()], depending depending on your desired
#' graph representation.
#'
#' ## Block memberships
#'
#' Let \eqn{Z_i} by a vector on the `k` dimensional simplex
#' representing the block memberships of node \eqn{i}.
#' To generate \eqn{z_i} we sample from a Dirichlet
#' distribution with parameter vector \eqn{\alpha}.
#' Block memberships for each node are independent.
#'
#' ## Degree heterogeneity
#'
#' In addition to block membership, the MMSBM also allows
#' nodes to have different propensities for edge formation.
#' We represent this propensity for node \eqn{i} by a positive
#' number \eqn{\theta_i}.
#'
#' ## Edge formulation
#'
#' Once we know the block membership vector \eqn{z_i, z_j} and the degree
#' heterogeneity parameters \eqn{\theta}, we need one more
#' ingredient, which is the baseline intensity of connections
#' between nodes in block `i` and block `j`. This is given by a
#' \eqn{k \times k} matrix \eqn{B}. Then each edge
#' \eqn{A_{i,j}} is Poisson distributed with parameter
#'
#' \deqn{
#'   \lambda_{i, j} = \theta_i \cdot z_i^T  B z_j \cdot \theta_j.
#' }
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_mmsbm <- mmsbm(n = 1000, k = 5, expected_density = 0.01)
#' lazy_mmsbm
#'
#' # sometimes you gotta let the world burn and
#' # sample a wildly dense graph
#'
#' dense_lazy_mmsbm <- mmsbm(n = 500, k = 3, expected_density = 0.8)
#' dense_lazy_mmsbm
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
#' alpha <- c(1, 2, 4, 1, 1)
#'
#' custom_mmsbm <- mmsbm(
#'   theta = theta,
#'   B = B,
#'   alpha = alpha,
#'   expected_degree = 50
#' )
#'
#' custom_mmsbm
#'
#' edgelist <- sample_edgelist(custom_mmsbm)
#' edgelist
#'
#' # efficient eigendecompostion that leverages low-rank structure in
#' # E(A) so that you don't have to form E(A) to find eigenvectors,
#' # as E(A) is typically dense. computation is
#' # handled via RSpectra
#'
#' population_eigs <- eigs_sym(custom_mmsbm)
#' svds(custom_mmsbm)$d
#'
mmsbm <- function(
    n = NULL, theta = NULL,
    k = NULL, B = NULL,
    ...,
    alpha = rep(1, k),
    sort_nodes = TRUE,
    force_pure = TRUE) {

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

  if (length(alpha) != nrow(B) || length(alpha) != ncol(B)) {
    stop("Both dimensions of B must match length of `alpha`", call. = FALSE)
  }

  # order mixing matrix by expected group size

  if (k > 1) {
    B <- B[order(alpha), ]
    B <- B[, order(alpha)]
  }

  alpha <- sort(alpha)

  Z <- t(igraph::sample_dirichlet(n, alpha))

  if (force_pure) {
    pure_indices <- sample(n, k)

    for (i in 1:k) {
      Z[pure_indices[i], ] <- 0
      Z[pure_indices[i], i] <- 1
    }
  }

  X <- Diagonal(n, theta) %*% Z

  if (sort_nodes) {
    idx <- sort_by_all_columns_idx(X)
    X <- X[idx, ]
    Z <- Z[idx, ]
    theta <- theta[idx]
  }

  mmsbm <- new_undirected_mmsbm(
    X = X,
    S = B,
    theta = theta,
    Z = Z,
    alpha = alpha,
    sorted = sort_nodes,
    ...
  )

  validate_undirected_mmsbm(mmsbm)
}

#' @method print undirected_mmsbm
#' @export
print.undirected_mmsbm <- function(x, ...) {

  cat(glue("Undirected Degree-Corrected Mixed Membership Stochastic Blockmodel\n", .trim = FALSE))
  cat(glue("------------------------------------------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k}\n\n", .trim = FALSE))

  cat("Traditional MMSBM parameterization:\n\n")
  cat("Block memberships portions (Z):", dim_and_class(x$Z), "\n")
  cat("Degree heterogeneity (theta):", dim_and_class(x$theta), "\n")
  cat("Block propensities (alpha):", dim_and_class(x$alpha), "\n\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
