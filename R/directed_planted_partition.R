validate_directed_planted_partition <- function(x) {
  values <- unclass(x)

  if (!inherits(x, "directed_planted_partition")) {
    stop(
      "Directed planted partition models must inherit ",
      "\"directed_planted_partition\" class!",
      call. = FALSE
    )
  }

  if (values$k1 != values$k2) {
    stop(
      "Directed planted partition models must have the same number of ",
      "outgoing and incoming blocks.",
      call. = FALSE
    )
  }

  distinct_theta_out_by_block <- vapply(
    split(values$theta_out, values$z_out),
    function(theta) length(unique(theta)),
    integer(1)
  )

  distinct_theta_in_by_block <- vapply(
    split(values$theta_in, values$z_in),
    function(theta) length(unique(theta)),
    integer(1)
  )

  if (any(distinct_theta_out_by_block > 1)) {
    stop(
      "`theta_out` must be equal for all nodes in a given block.",
      call. = FALSE
    )
  }

  if (any(distinct_theta_in_by_block > 1)) {
    stop(
      "`theta_in` must be equal for all nodes in a given block.",
      call. = FALSE
    )
  }

  x
}

check_scalar <- function(x, arg) {
  if (!is.null(x) && (!is.numeric(x) || length(x) != 1 || is.na(x))) {
    stop(glue("`{arg}` must be a scalar."), call. = FALSE)
  }

  invisible(x)
}

check_positive_integer <- function(x, arg) {
  if (
    !is.numeric(x) ||
      length(x) != 1 ||
      is.na(x) ||
      x < 1 ||
      x != round(x)
  ) {
    stop(glue("`{arg}` must be a positive integer."), call. = FALSE)
  }

  invisible(x)
}

#' Create a directed planted partition object
#'
#' To specify a directed planted partition model, you must specify
#' the number of nodes (via `n`), the number of blocks (via `k`), the
#' edge probabilities (either via `within_block` and `between_block`, or
#' via `a` and `b`), and the block sizes (optional, via `pi` or
#' `block_sizes`). By default, the communities a node sends edges from and
#' receives edges into are sampled independently; set `shared_blocks = TRUE`
#' to give each node a single community. We **strongly recommend**
#' setting the `expected_out_degree`, `expected_in_degree`, or
#' `expected_density` argument to avoid large memory allocations
#' associated with sampling large, dense graphs.
#'
#' @param n The number of nodes in the network. Must be
#'   a positive integer. This argument is required.
#'
#' @param k Number of planted partitions, as a positive integer. There are
#'   `k` outgoing blocks and `k` incoming blocks. Must be less than or equal
#'   to `n`. This argument is required.
#'
#' @param within_block Probability of a directed edge from a node in
#'   outgoing block `i` to a node in incoming block `i`, i.e. the diagonal
#'   of the mixing matrix `B`. Must be between zero and one. Must specify
#'   either `within_block` and `between_block`, or `a` and `b` to determine
#'   edge probabilities.
#'
#' @param between_block Probability of a directed edge from a node in
#'   outgoing block `i` to a node in incoming block `j`, for `i != j`, i.e.
#'   the off-diagonal of the mixing matrix `B`. Must be between zero and
#'   one. Must specify either `within_block` and `between_block`, or `a`
#'   and `b` to determine edge probabilities.
#'
#' @param a Number such that `a / n` is the probability of edges
#'   within a block. Useful for sparse graphs. Must specify either
#'   `within_block` and `between_block`, or `a` and `b` to determine
#'   edge probabilities.
#'
#' @param b Number such that `b / n` is the probability of edges
#'   between blocks. Useful for sparse graphs. Must specify either
#'   `within_block` and `between_block`, or `a` and `b` to determine
#'   edge probabilities.
#'
#' @param block_sizes (block sizes) Number of nodes in each block,
#'   as a vector of non-negative integers. Must have length `k` and
#'   must sum to `n`. Used for both the outgoing and the incoming blocks:
#'   exactly `block_sizes[j]` nodes are assigned to outgoing block `j`, and
#'   exactly `block_sizes[j]` nodes are assigned to incoming block `j`.
#'   Defaults to `NULL`. You can specify either `pi` or `block_sizes`,
#'   but not both.
#'
#' @param pi (block sizes) Relative block probabilities, used for both the
#'   outgoing and the incoming blocks. Must be positive, but do not need to
#'   sum to one, as they will be normalized internally. Must have length
#'   `k`. Defaults to `NULL`, in which case each node is equally likely to
#'   belong to each block (unless `block_sizes` is specified). You can
#'   specify either `pi` or `block_sizes`, but not both.
#'
#' @param shared_blocks Logical indicating whether each node should belong
#'   to the same outgoing and incoming block. Defaults to `FALSE`, in which
#'   case outgoing and incoming block memberships are sampled
#'   independently. When `TRUE`, a single block membership vector is
#'   sampled and used for both, so that `z_out` and `z_in` are identical.
#'   See [directed_dcsbm()] for details.
#'
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block. Useful for plotting.
#'   Defaults to `TRUE`. When `TRUE`, blocks are relabeled in increasing
#'   order of `pi` (or `block_sizes`), and the rows of the expected
#'   adjacency matrix are grouped by outgoing block while the columns are
#'   grouped by incoming block. Rows and columns are sorted independently
#'   when `shared_blocks = FALSE`, and together (so that row `i` and column
#'   `i` refer to the same node) when `shared_blocks = TRUE`.
#'
#' @inheritDotParams directed_factor_model expected_out_degree expected_in_degree expected_density
#' @inheritParams directed_factor_model
#'
#' @return A `directed_planted_partition` S3 object, which is a subclass
#'   of the [directed_dcsbm()] object (with `theta_out` and `theta_in`
#'   identically one), with additional fields:
#'
#'   - `within_block`: The probability of edge formation within a block,
#'     as specified (or as computed from `a / n`), before any rescaling
#'     by `expected_out_degree`, `expected_in_degree` or `expected_density`.
#'
#'   - `between_block`: The probability of edge formation between two
#'     distinct blocks, as specified (or as computed from `b / n`), before
#'     any rescaling by `expected_out_degree`, `expected_in_degree` or
#'     `expected_density`.
#'
#' @export
#' @family stochastic block models
#' @family directed graphs
#'
#' @details
#'
#' A directed planted partition model is a directed stochastic blockmodel
#' with `k` outgoing blocks and `k` incoming blocks, in which the diagonal
#' and the off-diagonal of the mixing matrix `B` are both constant. This
#' means that the probability of an edge from node \eqn{i} to node
#' \eqn{j} depends only on whether the outgoing block of \eqn{i} matches
#' the incoming block of \eqn{j}, but the particular blocks themselves
#' don't have any impact apart from this. It is a special case of
#' [directed_dcsbm()] where all degree-heterogeneity parameters are equal
#' to one.
#'
#' # Generative Model
#'
#' ## Block memberships
#'
#' Each node \eqn{i} has an outgoing block membership \eqn{z_{out, i}},
#' which governs the edges it sends, and an incoming block membership
#' \eqn{z_{in, i}}, which governs the edges it receives. Both take values
#' in `1, ..., k`. By default, outgoing memberships
#' are sampled independently across nodes from a categorical distribution
#' with parameter \eqn{\pi}, and incoming memberships are sampled in the
#' same way, independently of the outgoing memberships. When `block_sizes`
#' is specified, the memberships are instead a uniformly random assignment
#' of nodes to blocks with exactly `block_sizes[j]` nodes in block `j`
#' (again, independently for outgoing and incoming blocks).
#'
#' When `shared_blocks = TRUE`, a single membership vector is sampled
#' and \eqn{z_{out, i} = z_{in, i}} for every node.
#'
#' ## Mixing matrix
#'
#' The mixing matrix is the `k` by `k` matrix with `within_block` (or
#' `a / n`) on the diagonal and `between_block` (or `b / n`) everywhere
#' else, i.e.
#'
#' \deqn{
#'   B = p_{between} \mathbf{1} \mathbf{1}^T + (p_{within} - p_{between}) I_k,
#' }{
#'   B = p_between * 1 1^T + (p_within - p_between) * I_k,
#' }
#'
#' where \eqn{p_{within}} is `within_block` and \eqn{p_{between}} is
#' `between_block`. If you specify
#' `expected_out_degree`, `expected_in_degree` or `expected_density`, `B`
#' is rescaled by a constant, so only the ratio of `within_block` to
#' `between_block` matters.
#'
#' ## Edge formation
#'
#' Given the block memberships, each directed edge \eqn{A_{ij}} from node
#' \eqn{i} to node \eqn{j} is sampled independently (for every ordered
#' pair \eqn{(i, j)}, so the adjacency matrix is not symmetric) with
#' expectation
#'
#' \deqn{
#'   E[A_{ij}] = B_{z_{out, i}, z_{in, j}}.
#' }{
#'   E[A[i, j]] = B[z_out[i], z_in[j]].
#' }
#'
#' Edges are Poisson distributed when `poisson_edges = TRUE` (yielding a
#' multigraph) and Bernoulli distributed with the same probabilities when
#' `poisson_edges = FALSE`. Self-loops are allowed or not according to
#' `allow_self_loops`.
#'
#' ## Interpretation
#'
#' When `within_block > between_block`, nodes in outgoing block \eqn{s}
#' preferentially send edges to nodes in incoming block \eqn{s}. With the
#' default independent memberships, this is a co-clustering structure:
#' the communities that describe who a node sends edges to are unrelated
#' to the communities that describe who sends edges to it. With
#' `shared_blocks = TRUE`, each node has one community, and edges in both
#' directions are more common within communities than between them.
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_dpp <- directed_planted_partition(
#'   n = 100,
#'   k = 5,
#'   expected_density = 0.01,
#'   within_block = 0.1,
#'   between_block = 0.01
#' )
#'
#' lazy_dpp
#'
#' # sparse parameterization, with Bernoulli edges
#'
#' sparse_dpp <- directed_planted_partition(
#'   n = 1000,
#'   k = 3,
#'   a = 10,
#'   b = 4,
#'   poisson_edges = FALSE
#' )
#'
#' A <- sample_sparse(sparse_dpp)
#'
#' # only zeroes and ones!
#' range(A)
#'
#' # each node belongs to a single community, which determines both
#' # the edges it sends and the edges it receives
#'
#' shared_dpp <- directed_planted_partition(
#'   n = 100,
#'   k = 4,
#'   within_block = 0.3,
#'   between_block = 0.02,
#'   shared_blocks = TRUE
#' )
#'
#' identical(shared_dpp$z_out, shared_dpp$z_in)
#'
#' # exact block sizes
#'
#' sized_dpp <- directed_planted_partition(
#'   n = 90,
#'   k = 3,
#'   within_block = 0.3,
#'   between_block = 0.05,
#'   block_sizes = c(20, 30, 40)
#' )
#'
#' table(sized_dpp$z_out)
#' table(sized_dpp$z_in)
#'
#' population_svd <- svds(sized_dpp)
#' population_svd$d
#'
directed_planted_partition <- function(
  n,
  k,
  ...,
  within_block = NULL,
  between_block = NULL,
  a = NULL,
  b = NULL,
  block_sizes = NULL,
  pi = NULL,
  shared_blocks = FALSE,
  sort_nodes = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
) {
  check_positive_integer(n, "n")
  check_positive_integer(k, "k")

  if (k > n) {
    stop("`k` must be less than or equal to `n`.", call. = FALSE)
  }

  ### edge probabilities

  check_scalar(within_block, "within_block")
  check_scalar(between_block, "between_block")
  check_scalar(a, "a")
  check_scalar(b, "b")

  if (!is.null(within_block) && !is.null(a)) {
    stop("Must specify only one of `within_block` and `a`.", call. = FALSE)
  }

  if (!is.null(between_block) && !is.null(b)) {
    stop("Must specify only one of `between_block` and `b`.", call. = FALSE)
  }

  if (!is.null(a) && a < 0) {
    stop("`a` must be non-negative.", call. = FALSE)
  }

  if (!is.null(b) && b < 0) {
    stop("`b` must be non-negative.", call. = FALSE)
  }

  if (is.null(within_block) && !is.null(a)) {
    within_block <- a / n
  }

  if (is.null(between_block) && !is.null(b)) {
    between_block <- b / n
  }

  if (is.null(within_block) || is.null(between_block)) {
    stop(
      "Must specify either `within_block` and `between_block`, ",
      "or `a` and `b`.",
      call. = FALSE
    )
  }

  if (within_block < 0 || within_block > 1) {
    stop(
      "`within_block` (or `a / n`) must be between zero and one.",
      call. = FALSE
    )
  }

  if (between_block < 0 || between_block > 1) {
    stop(
      "`between_block` (or `b / n`) must be between zero and one.",
      call. = FALSE
    )
  }

  ### block sizes

  if (!is.null(block_sizes) && !is.null(pi)) {
    stop("Must specify at most one of `block_sizes` and `pi`.", call. = FALSE)
  }

  if (!is.null(block_sizes)) {
    check_block_sizes(
      block_sizes,
      arg = "block_sizes",
      k = k,
      n = n,
      k_description = "`k`",
      n_description = "`n`"
    )
  }

  if (!is.null(pi)) {
    if (length(pi) != k) {
      stop("Length of `pi` must match `k`.", call. = FALSE)
    }

    if (any(pi < 0)) {
      stop("All elements of `pi` must be >= 0.", call. = FALSE)
    }
  }

  B <- matrix(between_block, nrow = k, ncol = k)
  diag(B) <- within_block

  pp <- directed_dcsbm(
    theta_out = rep(1, n),
    theta_in = rep(1, n),
    B = B,
    ...,
    block_sizes_out = block_sizes,
    block_sizes_in = block_sizes,
    pi_out = pi,
    pi_in = pi,
    shared_blocks = shared_blocks,
    shared_theta = FALSE,
    sort_nodes = sort_nodes,
    force_identifiability = FALSE,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops,
    subclass = "directed_planted_partition"
  )

  pp$within_block <- within_block
  pp$between_block <- between_block

  validate_directed_planted_partition(pp)
}

#' @method print directed_planted_partition
#' @export
print.directed_planted_partition <- function(x, ...) {
  cat(glue("Directed Planted Partition Model\n", .trim = FALSE))
  cat(glue("--------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  shared_blocks <- if (isTRUE(x$shared_blocks)) {
    "TRUE (z_out and z_in are identical)"
  } else {
    "FALSE (z_out and z_in are independent)"
  }

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k1}\n", .trim = FALSE))
  cat(glue("Shared blocks: {shared_blocks}\n\n", .trim = FALSE))

  cat("Traditional planted partition parameterization:\n\n")
  cat("Within block edge probability:", signif(x$within_block, 4), "\n")
  cat("Between block edge probability:", signif(x$between_block, 4), "\n")
  cat("Block memberships (z_out):", dim_and_class(x$z_out), "\n")
  cat("Block memberships (z_in):", dim_and_class(x$z_in), "\n")
  cat("Block probabilities (pi_out):", dim_and_class(x$pi_out), "\n")
  cat("Block probabilities (pi_in):", dim_and_class(x$pi_in), "\n\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n")
  cat("Y:", dim_and_class(x$Y), "\n\n")

  cat("Poisson edges:", as.character(x$poisson_edges), "\n")
  cat("Allow self loops:", as.character(x$allow_self_loops), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue(
    "Expected out degree: {round(expected_out_degree(x), 1)}\n",
    .trim = FALSE
  ))
  cat(glue(
    "Expected in degree: {round(expected_in_degree(x), 1)}\n",
    .trim = FALSE
  ))
  cat(glue(
    "Expected density: {round(expected_density(x), 5)}\n",
    .trim = FALSE
  ))
}
