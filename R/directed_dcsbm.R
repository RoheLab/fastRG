new_directed_dcsbm <- function(
  X, S, Y,
  theta_in,
  theta_out,
  z_in,
  z_out,
  pi_in,
  pi_out,
  sorted,
  ...,
  subclass = character()) {
  subclass <- c(subclass, "directed_dcsbm")
  dcsbm <- directed_factor_model(X, S, Y, ..., subclass = subclass)
  dcsbm$theta_in <- theta_in
  dcsbm$theta_out <- theta_out
  dcsbm$z_in <- z_in
  dcsbm$z_out <- z_out
  dcsbm$pi_in <- pi_in
  dcsbm$pi_out <- pi_out
  dcsbm$sorted <- sorted
  dcsbm
}

validate_directed_dcsbm <- function(x) {

  values <- unclass(x)

  if (!is.factor(values$z_in)) {
    stop("`z_in` must be a factor.", call. = FALSE)
  }

  if (!is.factor(values$z_out)) {
    stop("`z_out` must be a factor.", call. = FALSE)
  }

  if (length(levels(values$z_in)) != values$k1) {
    stop(
      "Number of levels of `z1` must match the incoming rank of the model.",
      call. = FALSE
    )
  }

  if (length(levels(values$z_out)) != values$k2) {
    stop(
      "Number of levels of `z_out` must match the outgoing rank of the model.",
      call. = FALSE
    )
  }

  if (length(values$z_in) != nrow(values$X)) {
    stop("There must be one element of `z_in` for each row of `X`.")
  }

  if (length(values$z_out) != nrow(values$Y)) {
    stop("There must be one element of `z_out` for each row of `Y`.")
  }

  if (!is.numeric(values$theta_in)) {
    stop("`theta_in` must be a numeric.", call. = FALSE)
  }

  if (!is.numeric(values$theta_out)) {
    stop("`theta_out` must be a numeric.", call. = FALSE)
  }

  if (any(values$theta_in < 0)) {
    stop("Elements of `theta_in` must be strictly positive.", call. = FALSE)
  }

  if (any(values$theta_out < 0)) {
    stop("Elements of `theta_out` must be strictly positive.", call. = FALSE)
  }

  if (length(values$theta_in) != nrow(values$X)) {
    stop(
      "There must be one element of `theta_in` for each row of `X`.",
      call. = FALSE
    )
  }

  if (length(values$theta_out) != nrow(values$Y)) {
    stop(
      "There must be one element of `theta_out` for each row of `Y`.",
      call. = FALSE
    )
  }

  if (length(values$theta_in) != nrow(values$X)) {
    stop(
      "There must be one element of `theta_in` for each row of `X`.",
      call. = FALSE
    )
  }

  if (length(values$theta_out) != nrow(values$Y)) {
    stop(
      "There must be one element of `theta_out` for each row of `Y`.",
      call. = FALSE
    )
  }

  if (length(values$pi_in) != values$k1) {
    stop(
      "Dimension of `pi_in` must match the incoming rank of the model.",
      call. = FALSE
    )
  }

  if (length(values$pi_out) != values$k2) {
    stop(
      "Dimension of `pi_out` must match the outgoing rank of the model.",
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

#' Create a directed degree corrected stochastic blockmodel object
#'
#' To specify a degree-corrected stochastic blockmodel, you must specify
#' the degree-heterogeneity parameters (via `n_in` or `theta_in`, and
#' `n_out` or `theta_out`), the mixing matrix
#' (via `k_in` and `k_out`, or `B`), and the relative block
#' probabilities (optional, via `p_in` and `pi_out`).
#' We provide defaults for most of these
#' options to enable rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_in_degree`, `expected_out_degree`,
#' or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n (degree heterogeneity) The number of nodes in the blockmodel.
#'   Use when you don't want to specify the degree-heterogeneity
#'   parameters `theta_in` and `theta_out` by hand. When `n` is specified,
#'   `theta_in` and `theta_out` are randomly generated from
#'   a `LogNormal(2, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `n` defaults to `NULL`. You must specify either `n`
#'   or `theta_in` and `theta_out` together, but not both.
#'
#' @param theta_in (degree heterogeneity) A numeric vector
#'   explicitly specifying the degree heterogeneity
#'   parameters. This implicitly determines the number of nodes
#'   in the resulting graph, i.e. it will have `length(theta_in)` nodes.
#'   Must be positive. Setting to a vector of ones recovers
#'   a stochastic blockmodel without degree correction.
#'   Defaults to `NULL`. You must specify either `n`
#'   or `theta_in` and `theta_out` together, but not both.
#'
#' @param theta_out (degree heterogeneity) A numeric vector
#'   explicitly specifying the degree heterogeneity
#'   parameters. This implicitly determines the number of nodes
#'   in the resulting graph, i.e. it will have `length(theta)` nodes.
#'   Must be positive. Setting to a vector of ones recovers
#'   a stochastic blockmodel without degree correction.
#'   Defaults to `NULL`. You must specify either `n`
#'   or `theta_in` and `theta_out` together, but not both.
#'
#' @param k_in (mixing matrix) The number of blocks in the blockmodel.
#'   Use when you don't want to specify the mixing-matrix by hand.
#'   When `k_in` is specified, the elements of `B` are drawn
#'   randomly from a `Uniform(0, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `k_in` defaults to `NULL`. You must specify either `k_in` and
#'   `k_out` together, or `B`. You may specify all three at once, in which
#'   case `k_in` is only used to set `pi_in` (when `pi_in` is
#'   left at its default argument value).
#'
#' @param k_out (mixing matrix) The number of blocks in the blockmodel.
#'   Use when you don't want to specify the mixing-matrix by hand.
#'   When `k_out` is specified, the elements of `B` are drawn
#'   randomly from a `Uniform(0, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `k_out` defaults to `NULL`. You may specify all three at once, in which
#'   case `k_out` is only used to set `pi_out` (when `pi_out` is
#'   left at its default argument value).
#'
#' @param B (mixing matrix) A `k_in` by `k_out` matrix of block connection
#'   probabilities. The probability that a node in block `i` connects
#'   to a node in community `j` is `Poisson(B[i, j])`.
#'   `matrix` and `Matrix` objects are both acceptable.
#'   Defaults to `NULL`. You must specify either `k_in` and
#'   `k_out` together, or `B`, but not both.
#'
#' @param pi_in (relative block probabilities) Relative block
#'   probabilities. Must be positive, but do not need to sum
#'   to one, as they will be normalized internally.
#'   Must match the rows of `B`, or `k_in`. Defaults to
#'   `rep(1 / k_in, k_in)`, or a balanced incoming blocks.
#'
#' @param pi_out (relative block probabilities) Relative block
#'   probabilities. Must be positive, but do not need to sum
#'   to one, as they will be normalized internally.
#'   Must match the columns of `B`, or `k_out`. Defaults to
#'   `rep(1 / k_out, k_out)`, or a balanced outgoing blocks.
#'
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block. Useful for plotting.
#'   Defaults to `TRUE`.
#'
#' @param force_identifiability Logical indicating whether or not to
#'   normalize `theta_in` such that it sums to one within each incoming
#'   block and `theta_out` such that it sums to one within each outgoing
#'   block. Defaults to `TRUE`.
#'
#' @inheritDotParams directed_factor_model expected_in_degree expected_density expected_out_degree
#'
#' @return A `directed_dcsbm` S3 object, a subclass of the
#'   [directed_factor_model()] with the following additional
#'   fields:
#'
#'   - `theta_in`: A numeric vector of incoming community
#'     degree-heterogeneity parameters.
#'
#'   - `theta_out`: A numeric vector of outgoing community
#'     degree-heterogeneity parameters.
#'
#'   - `z_in`: The incoming community memberships of each node,
#'     as a [factor()]. The factor will have `k_in` levels,
#'     where `k_in` is the number of incoming
#'     communities in the stochastic blockmodel. There will not
#'     always necessarily be observed nodes in each community.
#'
#'   - `z_out`: The outgoing community memberships of each node,
#'     as a [factor()]. The factor will have `k_out` levels,
#'     where `k_out` is the number of outgoing
#'     communities in the stochastic blockmodel. There will not
#'     always necessarily be observed nodes in each community.
#'
#'   - `pi_in`: Sampling probabilities for each incoming
#'     community.
#'
#'   - `pi_out`: Sampling probabilities for each outgoing
#'     community.
#'
#'   - `sorted`: Logical indicating where nodes are arranged by
#'     block (and additionally by degree heterogeneity parameter)
#'     within each block.
#'
#' @export
#' @family stochastic block models
#' @family directed graphs
#'
#' @details
#'
#' # Generative Model
#'
#' There are two levels of randomness in a directed degree-corrected
#' stochastic blockmodel. First, we randomly chose a incoming
#' block membership and an outgoing block membership
#' for each node in the blockmodel. This is
#' handled by `directed_dcsbm()`. Then, given these block memberships,
#' we randomly sample edges between nodes. This second
#' operation is handled by [sample_edgelist()],
#' [sample_sparse()], [sample_igraph()] and
#' [sample_tidygraph()], depending on your desired
#' graph representation.
#'
#' ## Block memberships
#'
#' Let \eqn{x} represent the incoming block membership of a node
#' and \eqn{y} represent the outgoing block membership of a node.
#' To generate \eqn{x} we sample from a categorical
#' distribution with parameter \eqn{\pi_in}.
#' To generate \eqn{y} we sample from a categorical
#' distribution with parameter \eqn{\pi_out}.
#' Block memberships are independent across nodes. Incoming and outgoing
#' block memberships of the same node are also independent.
#'
#' ## Degree heterogeneity
#'
#' In addition to block membership, the DCSBM also
#' nodes to have different propensities for incoming and
#' outgoing edge formation.
#' We represent the propensity to form incoming edges for a
#' given node by a positive number \eqn{\theta_in}.
#' We represent the propensity to form outgoing edges for a
#' given node by a positive number \eqn{\theta_out}.
#' Typically the \eqn{\theta_in} (and \eqn{theta_out}) across all nodes are
#' constrained to sum to one for identifiability purposes,
#' but this doesn't really matter during sampling.
#'
#' ## Edge formulation
#'
#' Once we know the block memberships \eqn{x} and \eqn{y}
#' and the degree  heterogeneity parameters \eqn{\theta_in} and
#' \eqn{\theta_out}, we need one more
#' ingredient, which is the baseline intensity of connections
#' between nodes in block `i` and block `j`. Then each edge forms
#' independently according to a Poisson distribution with
#' parameters
#'
#' \deqn{
#'   \lambda = \theta_in * B[x, y] * \theta_out.
#' }
#'
#' @examples
#'
#' set.seed(27)
#'
#' B <- matrix(0.2, nrow = 5, ncol = 8)
#' diag(B) <- 0.9
#'
#' ddcsbm <- directed_dcsbm(
#'   n = 1000,
#'   B = B,
#'   k_in = 5,
#'   k_out = 8,
#'   expected_density = 0.01
#' )
#'
#' ddcsbm
#'
#' population_svd <- svds(ddcsbm)
#'
directed_dcsbm <- function(
  n = NULL,
  theta_in = NULL, theta_out = NULL,
  k_in = NULL, k_out = NULL, B = NULL,
  ...,
  pi_in = rep(1 / k_in, k_in),
  pi_out = rep(1 / k_out, k_out),
  sort_nodes = TRUE,
  force_identifiability = TRUE) {

  ### heterogeneity parameters

  if (is.null(n) && (is.null(theta_in) || is.null(theta_out))) {
    stop(
      "Must specify either `n` or both `theta_in` and `theta_out` together.",
      call. = FALSE
    )
  } else if (!is.null(n) && !is.null(theta_in) && !is.null(theta_out)) {
    stop(
      "Must specify either `n`, or both `theta_in` and `theta_out` together,",
      " but not all three at once.",
      call. = FALSE
    )
  } else if (is.null(theta_in) && is.null(theta_out)) {

    if (n < 1) {
      stop("`n` must be a positive integer.", call. = FALSE)
    }

    message(
      "Generating random degree heterogeneity parameters `theta_in` from a ",
      "LogNormal(2, 1) distribution. This distribution may change ",
      "in the future. Explicitly set `theta_in` for reproducible results.\n"
    )

    theta_in <- stats::rlnorm(n, meanlog = 2, sdlog = 1)
    theta_out <- stats::rlnorm(n, meanlog = 2, sdlog = 1)

  } else if (is.null(n)) {

    if (length(theta_in) != length(theta_out)) {
      stop(
        "Length of `theta_in` must match length of `theta_out`.",
        call. = FALSE
      )
    }

    n <- length(theta_in)
  }

  ### mixing matrix

  if (is.null(B) && (is.null(k_in) || is.null(k_out))) {
    stop(
      "Must specify either `k_in` and `k_out` together, or `B`.",
      call. = FALSE
    )
  } else if (is.null(B)) {

    if (k_in < 1) {
      stop("`k_in` must be a positive integer.", call. = FALSE)
    }

    if (k_out < 1) {
      stop("`k_out` must be a positive integer.", call. = FALSE)
    }

    message(
      "Generating random mixing matrix `B` with independent ",
      "Uniform(0, 1) entries. This distribution may change ",
      "in the future. Explicitly set `B` for reproducible results."
    )

    B <- Matrix(data = stats::runif(k_in * k_out), nrow = k_in, ncol = k_out)

  } else if (!is.null(B)) {

    k_in <- nrow(B)
    k_out <- ncol(B)
  }

  ### block membership

  if (length(pi_in) != nrow(B)) {
    stop("Length of `pi_in` must match number of rows in `B`.", call. = FALSE)
  }

  if (length(pi_out) != ncol(B)) {
    stop(
      "Length of `pi_out` must match number of columns in `B`.",
      call. = FALSE
    )
  }

  if (any(pi_in < 0)) {
    stop("All elements of `pi_in` must be >= 0.", call. = FALSE)
  }

  if (any(pi_out < 0)) {
    stop("All elements of `pi_out` must be >= 0.", call. = FALSE)
  }

  # order mixing matrix by expected group size

  if (k_in > 1) {
    B <- B[order(pi_in), ]
  }

  if (k_out > 1) {
    B <- B[, order(pi_out)]
  }

  pi_in <- sort(pi_in / sum(pi_in))
  pi_out <- sort(pi_out / sum(pi_out))

  # sample block memberships

  z_in <- sample(k_in, n, replace = TRUE, prob = pi_in)
  z_in <- factor(z_in, levels = 1:k_in, labels = paste0("incoming_block", 1:k_in))

  z_out <- sample(k_out, n, replace = TRUE, prob = pi_out)
  z_out <- factor(z_out, levels = 1:k_out, labels = paste0("outgoing_block", 1:k_out))

  if (sort_nodes) {
    z_in <- sort(z_in)
    z_out <- sort(z_out)
  }

  if (k_in > 1) {
    X <- sparse.model.matrix(~z_in + 0)
  } else {
    X <- Matrix(1, nrow = n, ncol = 1)
  }

  if (k_out > 1) {
    Y <- sparse.model.matrix(~z_out + 0)
  } else {
    Y <- Matrix(1, nrow = n, ncol = 1)
  }

  if (force_identifiability) {
    theta_in <- l1_normalize_within(theta_in, z_in)
    theta_out <- l1_normalize_within(theta_out, z_out)
  }


  X@x <- theta_in
  Y@x <- theta_out

  if (sort_nodes) {
    # X, Y indexing must match z indexing, be careful about this
    X <- sort_by_all_columns(X)
    Y <- sort_by_all_columns(Y)
  }

  if (k_in > 1) {
    theta_in <- X@x
  } else {
    theta_in <- sort(theta_in, decreasing = TRUE)
  }

  if (k_out > 1) {
    theta_out <- Y@x
  } else {
    theta_out <- sort(theta_out, decreasing = TRUE)
  }

  dcsbm <- new_directed_dcsbm(
    X = X,
    S = B,
    Y = Y,
    theta_in = theta_in,
    theta_out = theta_out,
    z_in = z_in,
    z_out = z_out,
    pi_in = pi_in,
    pi_out = pi_out,
    sorted = sort_nodes,
    ...
  )

  validate_directed_dcsbm(dcsbm)
}

#' @method print directed_dcsbm
#' @export
print.directed_dcsbm <- function(x, ...) {

  cat(glue("Directed Degree-Corrected Stochastic Blockmodel\n", .trim = FALSE))
  cat(glue("-----------------------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Incoming Blocks (k_in): {x$k1}\n", .trim = FALSE))
  cat(glue("Outgoing Blocks (k_out): {x$k2}\n\n", .trim = FALSE))

  cat("Traditional DCSBM parameterization:\n\n")
  cat("Block memberships (z_in):", dim_and_class(x$z_in), "\n")
  cat("Block memberships (z_out):", dim_and_class(x$z_out), "\n")
  cat("Degree heterogeneity (theta_in):", dim_and_class(x$theta_in), "\n")
  cat("Degree heterogeneity (theta_out):", dim_and_class(x$theta_out), "\n")
  cat("Block probabilities (pi_in):", dim_and_class(x$pi_in), "\n")
  cat("Block probabilities (pi_out):", dim_and_class(x$pi_out), "\n\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n")
  cat("Y:", dim_and_class(x$Y), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected in degree: {round(expected_in_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected out degree: {round(expected_out_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}\n", .trim = FALSE))
}
