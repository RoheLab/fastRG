new_directed_dcsbm <- function(
  X,
  S,
  Y,
  theta_out,
  theta_in,
  z_out,
  z_in,
  pi_out,
  pi_in,
  sorted,
  shared_blocks = FALSE,
  shared_theta = FALSE,
  ...,
  subclass = character()
) {
  subclass <- c(subclass, "directed_dcsbm")
  dcsbm <- directed_factor_model(X, S, Y, ..., subclass = subclass)
  dcsbm$theta_out <- theta_out
  dcsbm$theta_in <- theta_in
  dcsbm$z_out <- z_out
  dcsbm$z_in <- z_in
  dcsbm$pi_out <- pi_out
  dcsbm$pi_in <- pi_in
  dcsbm$sorted <- sorted
  dcsbm$shared_blocks <- shared_blocks
  dcsbm$shared_theta <- shared_theta
  dcsbm
}

validate_directed_dcsbm <- function(x) {
  values <- unclass(x)

  if (!is.factor(values$z_out)) {
    stop("`z_out` must be a factor.", call. = FALSE)
  }

  if (!is.factor(values$z_in)) {
    stop("`z_in` must be a factor.", call. = FALSE)
  }

  if (length(levels(values$z_out)) != values$k1) {
    stop(
      "Number of levels of `z1` must match the incoming rank of the model.",
      call. = FALSE
    )
  }

  if (length(levels(values$z_in)) != values$k2) {
    stop(
      "Number of levels of `z_in` must match the outgoing rank of the model.",
      call. = FALSE
    )
  }

  if (length(values$z_out) != nrow(values$X)) {
    stop("There must be one element of `z_out` for each row of `X`.")
  }

  if (length(values$z_in) != nrow(values$Y)) {
    stop("There must be one element of `z_in` for each row of `Y`.")
  }

  if (!is.numeric(values$theta_out)) {
    stop("`theta_out` must be a numeric.", call. = FALSE)
  }

  if (!is.numeric(values$theta_in)) {
    stop("`theta_in` must be a numeric.", call. = FALSE)
  }

  if (any(values$theta_out < 0)) {
    stop("Elements of `theta_out` must be strictly positive.", call. = FALSE)
  }

  if (any(values$theta_in < 0)) {
    stop("Elements of `theta_in` must be strictly positive.", call. = FALSE)
  }

  if (length(values$theta_out) != nrow(values$X)) {
    stop(
      "There must be one element of `theta_out` for each row of `X`.",
      call. = FALSE
    )
  }

  if (length(values$theta_in) != nrow(values$Y)) {
    stop(
      "There must be one element of `theta_in` for each row of `Y`.",
      call. = FALSE
    )
  }

  if (length(values$theta_out) != nrow(values$X)) {
    stop(
      "There must be one element of `theta_out` for each row of `X`.",
      call. = FALSE
    )
  }

  if (length(values$theta_in) != nrow(values$Y)) {
    stop(
      "There must be one element of `theta_in` for each row of `Y`.",
      call. = FALSE
    )
  }

  if (length(values$pi_out) != values$k1) {
    stop(
      "Dimension of `pi_out` must match the incoming rank of the model.",
      call. = FALSE
    )
  }

  if (length(values$pi_in) != values$k2) {
    stop(
      "Dimension of `pi_in` must match the outgoing rank of the model.",
      call. = FALSE
    )
  }

  if (min(values$S) < 0) {
    stop(
      "All elements of `B` must be non-negative.",
      call. = FALSE
    )
  }

  if (isTRUE(values$shared_blocks) && !identical(values$z_out, values$z_in)) {
    stop(
      "`z_out` and `z_in` must be identical when `shared_blocks = TRUE`.",
      call. = FALSE
    )
  }

  x
}

# check that `block_sizes` is a valid vector of block counts for `k` blocks
# and `n` nodes. `arg`, `k_description` and `n_description` only appear in
# error messages
check_block_sizes <- function(
  block_sizes,
  arg,
  k,
  n,
  k_description,
  n_description
) {
  if (
    !is.numeric(block_sizes) ||
      anyNA(block_sizes) ||
      any(block_sizes < 0) ||
      any(block_sizes != round(block_sizes))
  ) {
    stop(
      glue("All elements of `{arg}` must be non-negative integers."),
      call. = FALSE
    )
  }

  if (length(block_sizes) != k) {
    stop(
      glue("Length of `{arg}` must match {k_description}."),
      call. = FALSE
    )
  }

  if (sum(block_sizes) != n) {
    stop(
      glue("Sum of `{arg}` must equal {n_description}."),
      call. = FALSE
    )
  }

  invisible(block_sizes)
}

# sample block memberships for `n` nodes in `k` blocks. when `block_sizes` is
# `NULL`, memberships are drawn independently with probabilities `pi`.
# otherwise exactly `block_sizes[j]` nodes are assigned to block `j`, and
# which nodes land in which block is uniformly random
sample_block_memberships <- function(n, k, pi, block_sizes = NULL) {
  if (is.null(block_sizes)) {
    return(sample(k, n, replace = TRUE, prob = pi))
  }

  # index rather than calling sample() on the vector directly, since
  # sample(x) samples from 1:x when length(x) == 1
  memberships <- rep(seq_len(k), times = block_sizes)
  memberships[sample.int(length(memberships))]
}

#' Create a directed degree corrected stochastic blockmodel object
#'
#' To specify a degree-corrected stochastic blockmodel, you must specify
#' the degree-heterogeneity parameters (via `n` or `theta_out` and
#' `theta_in`), the mixing matrix
#' (via `k_out` and `k_in`, or `B`), and the block sizes (optional,
#' either as relative block probabilities via `pi_out` and `pi_in`, or as
#' exact block counts via `block_sizes_out` and `block_sizes_in`).
#' By default, outgoing and incoming block memberships and degree
#' parameters are unrelated (a co-factor model); use `shared_blocks` and
#' `shared_theta` to tie them together.
#' We provide defaults for most of these
#' options to enable rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_out_degree`, `expected_in_degree`,
#' or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n (degree heterogeneity) The number of nodes in the blockmodel.
#'   Use when you don't want to specify the degree-heterogeneity
#'   parameters `theta_out` and `theta_in` by hand. When `n` is specified,
#'   `theta_out` and `theta_in` are randomly generated from
#'   a `LogNormal(2, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `n` defaults to `NULL`. You must specify either `n`
#'   or `theta_out` and `theta_in` together, but not both.
#'   When `shared_theta = TRUE`, a single vector is generated and used
#'   for both `theta_out` and `theta_in`.
#'
#' @param theta_out (degree heterogeneity) A numeric vector
#'   explicitly specifying the degree heterogeneity
#'   parameters. This implicitly determines the number of nodes
#'   in the resulting graph, i.e. it will have `length(theta_out)` nodes.
#'   Must be positive. Setting to a vector of ones recovers
#'   a stochastic blockmodel without degree correction.
#'   Defaults to `NULL`. You must specify either `n`
#'   or `theta_out` and `theta_in` together, but not both. `theta_out`
#'   controls outgoing degree propensity, or, equivalently,
#'   row sums of the adjacency matrix.
#'
#' @param theta_in (degree heterogeneity) A numeric vector
#'   explicitly specifying the degree heterogeneity
#'   parameters. This implicitly determines the number of nodes
#'   in the resulting graph, i.e. it will have `length(theta)` nodes.
#'   Must be positive. Setting to a vector of ones recovers
#'   a stochastic blockmodel without degree correction.
#'   Defaults to `NULL`. You must specify either `n`
#'   or `theta_out` and `theta_in` together, but not both. `theta_in`
#'   controls incoming degree propensity, or, equivalently, column
#'   sums of the adjacency matrix.
#'
#' @param k_out (mixing matrix) The number of outgoing blocks in the blockmodel.
#'   Use when you don't want to specify the mixing-matrix by hand.
#'   When `k_out` is specified, the elements of `B` are drawn
#'   randomly from a `Uniform(0, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `k_out` defaults to `NULL`. You must specify either `k_out` and
#'   `k_in` together, or `B`. You may specify all three at once, in which
#'   case the number of outgoing blocks is taken from `B`.
#'
#' @param k_in (mixing matrix) The number of incoming blocks in the blockmodel.
#'   Use when you don't want to specify the mixing-matrix by hand.
#'   When `k_in` is specified, the elements of `B` are drawn
#'   randomly from a `Uniform(0, 1)` distribution.
#'   This is subject to change, and may not be reproducible.
#'   `k_in` defaults to `NULL`. You may specify all three at once, in which
#'   case the number of incoming blocks is taken from `B`.
#'
#' @param B (mixing matrix) A `k_out` by `k_in` matrix of block connection
#'   probabilities. The probability that a node in block `i` connects
#'   to a node in community `j` is `Poisson(B[i, j])`.
#'   `matrix` and `Matrix` objects are both acceptable.
#'   Defaults to `NULL`. You must specify either `k_out` and
#'   `k_in` together, or `B`, but not both.
#'
#' @param block_sizes_out (block sizes) Number of nodes in each outgoing
#'   block, as a vector of non-negative integers. Must match the number of
#'   rows of `B` (or `k_out`) and must sum to `n` (or `length(theta_out)`).
#'   When specified, exactly `block_sizes_out[j]` nodes are assigned to
#'   outgoing block `j`, and which nodes land in which block is uniformly
#'   random. Defaults to `NULL`, in which case outgoing blocks are made to
#'   be as balanced as possible, unless `pi_out` is specified. You can
#'   specify either `pi_out` or `block_sizes_out`, but not both.
#'
#' @param block_sizes_in (block sizes) Number of nodes in each incoming
#'   block, as a vector of non-negative integers. Must match the number of
#'   columns of `B` (or `k_in`) and must sum to `n` (or `length(theta_in)`).
#'   When specified, exactly `block_sizes_in[j]` nodes are assigned to
#'   incoming block `j`, and which nodes land in which block is uniformly
#'   random. Defaults to `NULL`, in which case incoming blocks are made to
#'   be as balanced as possible, unless `pi_in` is specified. You can
#'   specify either `pi_in` or `block_sizes_in`, but not both.
#'
#' @param pi_out (block sizes) Relative block
#'   probabilities for outgoing blocks. Must be positive, but do not need to
#'   sum to one, as they will be normalized internally.
#'   Must match the rows of `B`, or `k_out`. When specified, outgoing block
#'   memberships are sampled independently for each node with these
#'   probabilities. Defaults to `NULL`, in which case the `block_sizes_out`
#'   argument (or its balanced default) takes precedence. You can specify
#'   either `pi_out` or `block_sizes_out`, but not both.
#'
#' @param pi_in (block sizes) Relative block
#'   probabilities for incoming blocks. Must be positive, but do not need to
#'   sum to one, as they will be normalized internally.
#'   Must match the columns of `B`, or `k_in`. When specified, incoming
#'   block memberships are sampled independently for each node with these
#'   probabilities. Defaults to `NULL`, in which case the `block_sizes_in`
#'   argument (or its balanced default) takes precedence. You can specify
#'   either `pi_in` or `block_sizes_in`, but not both.
#'
#' @param shared_blocks (block memberships) Logical indicating whether each
#'   node should belong to the same outgoing and incoming block. Defaults to
#'   `FALSE`, in which case outgoing and incoming block memberships are
#'   sampled independently. When `TRUE`, a single block membership vector is
#'   sampled and used for both, so that `z_out` and `z_in` are identical.
#'   This requires the same number of outgoing and incoming blocks, i.e.
#'   `B` must be square (or `k_out` must equal `k_in`; specifying only one
#'   of them is enough). It is enough to specify only one of `pi_out` and
#'   `pi_in` (or of `block_sizes_out` and `block_sizes_in`); if both are
#'   specified, they must agree.
#'
#' @param shared_theta (degree heterogeneity) Logical indicating whether
#'   each node should have the same propensity to form outgoing and incoming
#'   edges. Defaults to `FALSE`, in which case `theta_out` and `theta_in`
#'   are unrelated. When `TRUE`, a single degree-heterogeneity vector is
#'   used for both `theta_out` and `theta_in`: either specify `n` (one
#'   vector is then generated at random), or only one of `theta_out` and
#'   `theta_in`. If you specify both `theta_out` and `theta_in`, they must
#'   be identical. Requires `shared_blocks = TRUE`, so that the stored
#'   `theta_out` and `theta_in` are identical (also after normalization
#'   by `force_identifiability`).
#'
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block. Useful for plotting.
#'   Defaults to `TRUE`. When `TRUE`, rows of the expected adjacency matrix
#'   are first sorted by outgoing block membership, and then by
#'   outgoing degree-correction parameters within each outgoing block.
#'   A similar sorting procedure occurs independently for the columns,
#'   according to the incoming blocks.
#'   Additionally, `pi_out` and `pi_in` are sorted in increasing order,
#'   and the rows and columns of the `B` matrix are permuted to match the
#'   new orderings. When `block_sizes_out` or `block_sizes_in` are
#'   specified, they are permuted in the same way, so each block keeps
#'   exactly the requested number of nodes and its row (or column) of `B`.
#'
#'   **Rows and columns index different nodes.** In the default co-factor
#'   model (`shared_blocks = FALSE`), the rows and the columns are permuted
#'   *independently*, so after sorting, row `i` and column `i` of the
#'   adjacency matrix are generally not the same node. In particular, the
#'   diagonal of `A` (self-loops, see `allow_self_loops`) does not mean
#'   "node `i` to itself", and `X[i, ]` and `Y[i, ]` describe different
#'   nodes. If you need rows and columns to index the same nodes, use
#'   `shared_blocks = TRUE` or `sort_nodes = FALSE`.
#'
#'   When `shared_blocks = TRUE`, rows and columns are instead permuted
#'   together, so that row `i` and column `i` of the expected adjacency
#'   matrix refer to the same node: nodes are sorted by block, and then by
#'   decreasing `theta_out` and `theta_in` within each block.
#'
#' @param force_identifiability Logical indicating whether or not to
#'   normalize `theta_out` such that it sums to one within each incoming
#'   block and `theta_in` such that it sums to one within each outgoing
#'   block. Defaults to `TRUE`.
#'
#' @inheritDotParams directed_factor_model expected_out_degree expected_density expected_in_degree
#' @inheritParams directed_factor_model
#'
#' @return A `directed_dcsbm` S3 object, a subclass of the
#'   [directed_factor_model()] with the following additional
#'   fields:
#'
#'   - `theta_out`: A numeric vector of outgoing
#'     degree-heterogeneity parameters.
#'
#'   - `theta_in`: A numeric vector of incoming
#'     degree-heterogeneity parameters.
#'
#'   - `z_out`: The outgoing community memberships of each node,
#'     as a [factor()]. The factor will have `k_out` levels,
#'     where `k_out` is the number of outgoing
#'     communities in the stochastic blockmodel. There will not
#'     always necessarily be observed nodes in each community.
#'
#'   - `z_in`: The incoming community memberships of each node,
#'     as a [factor()]. The factor will have `k_in` levels,
#'     where `k_in` is the number of incoming
#'     communities in the stochastic blockmodel. There will not
#'     always necessarily be observed nodes in each community.
#'
#'   - `pi_out`: Sampling probabilities for each outgoing
#'     community.
#'
#'   - `pi_in`: Sampling probabilities for each incoming
#'     community.
#'
#'   - `sorted`: Logical indicating where nodes are arranged by
#'     block (and additionally by degree heterogeneity parameter)
#'     within each block.
#'
#'   - `shared_blocks`: Logical indicating whether outgoing and incoming
#'     block memberships are shared, in which case `z_out` and `z_in`
#'     are identical.
#'
#'   - `shared_theta`: Logical indicating whether outgoing and incoming
#'     degree-heterogeneity parameters are shared, in which case `theta_out`
#'     and `theta_in` are identical.
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
#' distribution with parameter \eqn{\pi_{out}}.
#' To generate \eqn{y} we sample from a categorical
#' distribution with parameter \eqn{\pi_{in}}.
#' Block memberships are independent across nodes. Incoming and outgoing
#' block memberships of the same node are also independent.
#'
#' These categorical draws are used when `pi_out` (`pi_in`) is specified.
#' Otherwise, the outgoing (incoming) block memberships are a uniformly
#' random assignment of nodes to blocks with fixed block sizes: either
#' exactly the sizes given in `block_sizes_out` (`block_sizes_in`), or, by
#' default, block sizes that are as balanced as possible (every block gets
#' `floor(n / k_out)` nodes, and the remainder is spread over the first
#' blocks, one node each), as in [dcsbm()]. Fixed block sizes are useful,
#' for instance, to construct blockmodels with exactly repeated singular
#' values.
#'
#' ## Degree heterogeneity
#'
#' In addition to block membership, the DCSBM also
#' nodes to have different propensities for incoming and
#' outgoing edge formation.
#' We represent the propensity to form incoming edges for a
#' given node by a positive number \eqn{\theta_{out}}.
#' We represent the propensity to form outgoing edges for a
#' given node by a positive number \eqn{\theta_{in}}.
#' Typically the \eqn{\theta_{out}} (and \eqn{\theta_{in}}) across all nodes are
#' constrained to sum to one for identifiability purposes,
#' but this doesn't really matter during sampling.
#'
#' ## Edge formulation
#'
#' Once we know the block memberships \eqn{x} and \eqn{y}
#' and the degree  heterogeneity parameters \eqn{\theta_{in}} and
#' \eqn{\theta_{out}}, we need one more
#' ingredient, which is the baseline intensity of connections
#' between nodes in block `i` and block `j`. Then each edge forms
#' independently according to a Poisson distribution with
#' parameters
#'
#' \deqn{
#'   \lambda = \theta_{in} \cdot B_{x, y} \cdot \theta_{out}.
#' }
#'
#' # Co-factor models versus shared blocks and degree heterogeneity
#'
#' By default (`shared_blocks = FALSE` and `shared_theta = FALSE`),
#' `directed_dcsbm()` is a *co-factor* model: every node has an outgoing
#' block \eqn{z_{out}} and degree parameter \eqn{\theta_{out}} that govern
#' the edges it sends (the rows of \eqn{A}), and a separate incoming block
#' \eqn{z_{in}} and degree parameter \eqn{\theta_{in}} that govern the edges
#' it receives (the columns of \eqn{A}). These are drawn independently, so
#' the communities a node sends edges from and receives edges into are
#' unrelated, and the two sets of communities can even have different
#' numbers of blocks. This is the model targeted by co-clustering methods
#' that estimate separate sending and receiving communities.
#'
#' In the co-factor model, `sort_nodes = TRUE` (the default) sorts the rows
#' and the columns of \eqn{A} *independently*, so row \eqn{i} and column
#' \eqn{i} are generally not the same node. Consequently the diagonal of
#' \eqn{A}, and so the self-loops controlled by `allow_self_loops`, is not
#' meaningful as "node \eqn{i} to itself": entry \eqn{A_{ii}} connects the
#' \eqn{i}-th sender to the \eqn{i}-th receiver in their respective
#' orderings. If you need rows and columns to index the same nodes (for
#' example, to interpret self-loops or to read the graph as a directed
#' graph on a single set of nodes), use `shared_blocks = TRUE` or
#' `sort_nodes = FALSE`.
#'
#' Often it is more natural to think of each node as belonging to a single
#' community, which shapes both the edges it sends and the edges it
#' receives. Setting `shared_blocks = TRUE` draws one membership vector
#' \eqn{z} with \eqn{z_{out} = z_{in} = z}, so that
#' \eqn{E[A_{ij}] = \theta_{out, i} B_{z_i, z_j} \theta_{in, j}}, where
#' an asymmetric `B` encodes directed relationships between the communities.
#' With shared blocks you can additionally set `shared_theta = TRUE`, which
#' uses one degree-heterogeneity vector \eqn{\theta} for both
#' \eqn{\theta_{out}} and \eqn{\theta_{in}}, so nodes that send many edges
#' also tend to receive many edges. (`shared_theta = TRUE` requires
#' `shared_blocks = TRUE`.) With shared blocks, `sort_nodes = TRUE`
#' permutes rows and columns together, so that row `i` and column `i` of
#' \eqn{A} always refer to the same node.
#'
#' @examples
#'
#' set.seed(27)
#'
#' B <- matrix(0.2, nrow = 5, ncol = 8)
#' diag(B) <- 0.9
#'
#' ddcsbm <- directed_dcsbm(
#'   n = 100,
#'   B = B,
#'   k_out = 5,
#'   k_in = 8,
#'   expected_density = 0.01
#' )
#'
#' ddcsbm
#'
#' population_svd <- svds(ddcsbm)
#'
#' # specify the exact number of nodes in each outgoing and incoming block,
#' # rather than sampling block memberships
#'
#' ddcsbm_explicit_block_sizes <- directed_dcsbm(
#'   theta_out = rexp(100, 1 / 3) + 1,
#'   theta_in = rexp(100, 1 / 3) + 1,
#'   B = B,
#'   block_sizes_out = c(10, 15, 20, 25, 30),
#'   block_sizes_in = c(5, 5, 10, 10, 15, 15, 20, 20),
#'   expected_density = 0.01
#' )
#'
#' # respects block sizes
#' table(ddcsbm_explicit_block_sizes$z_out)
#' table(ddcsbm_explicit_block_sizes$z_in)
#'
#' # each node belongs to the same outgoing and incoming community, and has
#' # the same propensity to send and receive edges
#'
#' B_square <- matrix(0.1, nrow = 4, ncol = 4)
#' diag(B_square) <- 0.8
#' B_square[1, 2] <- 0.5 # block 1 sends many edges to block 2
#'
#' shared_ddcsbm <- directed_dcsbm(
#'   theta_out = rexp(100, 1 / 3) + 1,
#'   B = B_square,
#'   shared_blocks = TRUE,
#'   shared_theta = TRUE,
#'   expected_density = 0.05
#' )
#'
#' shared_ddcsbm
#'
#' identical(shared_ddcsbm$z_out, shared_ddcsbm$z_in)
#'
directed_dcsbm <- function(
  n = NULL,
  theta_out = NULL,
  theta_in = NULL,
  k_out = NULL,
  k_in = NULL,
  B = NULL,
  ...,
  block_sizes_out = NULL,
  block_sizes_in = NULL,
  pi_out = NULL,
  pi_in = NULL,
  shared_blocks = FALSE,
  shared_theta = FALSE,
  sort_nodes = TRUE,
  force_identifiability = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
) {
  # NOTE:
  #   - X corresponds to outgoing communities, outgoing edges and rows of A
  #   - Y corresponds to incoming communities, incoming edges and columns of A

  check_flag(shared_blocks, "shared_blocks")
  check_flag(shared_theta, "shared_theta")

  if (shared_theta && !shared_blocks) {
    stop(
      "`shared_theta = TRUE` requires `shared_blocks = TRUE`. Shared degree ",
      "heterogeneity is only supported when outgoing and incoming blocks ",
      "are shared as well.",
      call. = FALSE
    )
  }

  ### heterogeneity parameters

  if (shared_theta) {
    if (
      !is.null(theta_out) &&
        !is.null(theta_in) &&
        (length(theta_out) != length(theta_in) || any(theta_out != theta_in))
    ) {
      stop(
        "When `shared_theta = TRUE`, `theta_out` and `theta_in` must be ",
        "identical if both are specified. Specify only one of them.",
        call. = FALSE
      )
    }

    theta <- if (is.null(theta_out)) theta_in else theta_out

    if (is.null(n) && is.null(theta)) {
      stop(
        "Must specify either `n`, or `theta_out` or `theta_in`.",
        call. = FALSE
      )
    } else if (!is.null(n) && !is.null(theta)) {
      stop(
        "Must specify only one of `n` and `theta_out` / `theta_in`, ",
        "not both.",
        call. = FALSE
      )
    } else if (is.null(theta)) {
      if (n < 1) {
        stop("`n` must be a positive integer.", call. = FALSE)
      }

      message(
        "Generating random degree heterogeneity parameters `theta_out` and `theta_in` from ",
        "a single LogNormal(2, 1) draw, shared by `theta_out` and `theta_in`. ",
        "This distribution may change in the future. Explicitly set `theta_out` ",
        "or `theta_in` for reproducible results.\n"
      )

      theta <- stats::rlnorm(n, meanlog = 2, sdlog = 1)
    } else {
      n <- length(theta)
    }

    theta_out <- theta
    theta_in <- theta
  } else if (is.null(n) && (is.null(theta_out) || is.null(theta_in))) {
    stop(
      "Must specify either `n` or both `theta_out` and `theta_in` together.",
      call. = FALSE
    )
  } else if (!is.null(n) && !is.null(theta_out) && !is.null(theta_in)) {
    stop(
      "Must specify either `n`, or both `theta_out` and `theta_in` together,",
      " but not all three at once.",
      call. = FALSE
    )
  } else if (is.null(theta_out) && is.null(theta_in)) {
    if (n < 1) {
      stop("`n` must be a positive integer.", call. = FALSE)
    }

    message(
      "Generating random degree heterogeneity parameters `theta_out` and `theta_in` from ",
      "LogNormal(2, 1) distributions. This distribution may change ",
      "in the future. Explicitly set `theta_out` and `theta_in` for reproducible results.\n"
    )

    theta_out <- stats::rlnorm(n, meanlog = 2, sdlog = 1)
    theta_in <- stats::rlnorm(n, meanlog = 2, sdlog = 1)
  } else if (is.null(n)) {
    if (length(theta_out) != length(theta_in)) {
      stop(
        "Length of `theta_out` must match length of `theta_in`.",
        call. = FALSE
      )
    }

    n <- length(theta_out)
  } else {
    stop(
      "Must specify either `n`, or both `theta_out` and `theta_in` together,",
      " but not both.",
      call. = FALSE
    )
  }

  ### mixing matrix

  if (shared_blocks && is.null(B)) {
    # with shared blocks there is only one number of blocks, so it is
    # enough to specify either `k_out` or `k_in`
    if (is.null(k_out)) k_out <- k_in
    if (is.null(k_in)) k_in <- k_out

    if (!is.null(k_out) && k_out != k_in) {
      stop(
        "`k_out` and `k_in` must be equal when `shared_blocks = TRUE`.",
        call. = FALSE
      )
    }
  }

  if (shared_blocks && !is.null(B) && nrow(B) != ncol(B)) {
    stop(
      "`B` must be a square matrix when `shared_blocks = TRUE`.",
      call. = FALSE
    )
  }

  if (is.null(B) && (is.null(k_out) || is.null(k_in))) {
    stop(
      "Must specify either `k_out` and `k_in` together, or `B`.",
      call. = FALSE
    )
  } else if (is.null(B)) {
    if (k_out < 1) {
      stop("`k_out` must be a positive integer.", call. = FALSE)
    }

    if (k_in < 1) {
      stop("`k_in` must be a positive integer.", call. = FALSE)
    }

    message(
      "Generating random mixing matrix `B` with independent ",
      "Uniform(0, 1) entries. This distribution may change ",
      "in the future. Explicitly set `B` for reproducible results."
    )

    B <- Matrix(data = stats::runif(k_out * k_in), nrow = k_out, ncol = k_in)
  } else if (!is.null(B)) {
    k_out <- nrow(B)
    k_in <- ncol(B)
  }

  ### block membership

  if (shared_blocks) {
    # a single membership vector is used for both outgoing and incoming
    # blocks, so the outgoing and incoming specifications must agree. it is
    # enough to specify only one of them
    if (
      !is.null(block_sizes_out) &&
        !is.null(block_sizes_in) &&
        (length(block_sizes_out) != length(block_sizes_in) ||
          any(block_sizes_out != block_sizes_in))
    ) {
      stop(
        "When `shared_blocks = TRUE`, `block_sizes_out` and `block_sizes_in` ",
        "must be identical if both are specified. Specify only one of them.",
        call. = FALSE
      )
    }

    if (
      !is.null(pi_out) &&
        !is.null(pi_in) &&
        (length(pi_out) != length(pi_in) ||
          !isTRUE(all.equal(pi_out / sum(pi_out), pi_in / sum(pi_in))))
    ) {
      stop(
        "When `shared_blocks = TRUE`, `pi_out` and `pi_in` must agree ",
        "(up to normalization) if both are specified. Specify only one of them.",
        call. = FALSE
      )
    }

    block_sizes <- if (is.null(block_sizes_out)) block_sizes_in else block_sizes_out
    pi <- if (is.null(pi_out)) pi_in else pi_out

    if (!is.null(block_sizes) && !is.null(pi)) {
      stop(
        "When `shared_blocks = TRUE`, must specify at most one of block sizes ",
        "(`block_sizes_out` / `block_sizes_in`) and block probabilities ",
        "(`pi_out` / `pi_in`).",
        call. = FALSE
      )
    }

    block_sizes_out <- block_sizes
    block_sizes_in <- block_sizes
    pi_out <- pi
    pi_in <- pi
  }

  if (!is.null(block_sizes_out) && !is.null(pi_out)) {
    stop(
      "Must specify at most one of `block_sizes_out` and `pi_out`.",
      call. = FALSE
    )
  }

  if (!is.null(block_sizes_in) && !is.null(pi_in)) {
    stop(
      "Must specify at most one of `block_sizes_in` and `pi_in`.",
      call. = FALSE
    )
  }

  if (!is.null(block_sizes_out)) {
    check_block_sizes(
      block_sizes_out,
      arg = "block_sizes_out",
      k = nrow(B),
      n = n,
      k_description = "number of rows in `B`",
      n_description = "`n` or `length(theta_out)`"
    )

    # for sorting by block size later
    pi_out <- block_sizes_out / n
  } else if (is.null(pi_out)) {
    # default: blocks as balanced as possible, like `dcsbm()`
    block_sizes_out <- balanced_block_sizes(n, k_out)
    pi_out <- block_sizes_out / n
  }

  if (!is.null(block_sizes_in)) {
    check_block_sizes(
      block_sizes_in,
      arg = "block_sizes_in",
      k = ncol(B),
      n = n,
      k_description = "number of columns in `B`",
      n_description = "`n` or `length(theta_in)`"
    )

    # for sorting by block size later
    pi_in <- block_sizes_in / n
  } else if (is.null(pi_in)) {
    # default: blocks as balanced as possible, like `dcsbm()`
    block_sizes_in <- balanced_block_sizes(n, k_in)
    pi_in <- block_sizes_in / n
  }

  if (length(pi_out) != nrow(B)) {
    stop("Length of `pi_out` must match number of rows in `B`.", call. = FALSE)
  }

  if (length(pi_in) != ncol(B)) {
    stop(
      "Length of `pi_in` must match number of columns in `B`.",
      call. = FALSE
    )
  }

  if (any(pi_out < 0)) {
    stop("All elements of `pi_out` must be >= 0.", call. = FALSE)
  }

  if (any(pi_in < 0)) {
    stop("All elements of `pi_in` must be >= 0.", call. = FALSE)
  }

  # order mixing matrix by expected group size. when block sizes are
  # specified they are permuted alongside, so that block labels, the
  # rows / columns of B, pi and the block sizes all stay consistent

  if (k_out > 1 && sort_nodes) {
    out_order <- order(pi_out)
    B <- B[out_order, , drop = FALSE]
    pi_out <- pi_out[out_order]
    block_sizes_out <- block_sizes_out[out_order]
  }

  if (k_in > 1 && sort_nodes) {
    in_order <- order(pi_in)
    B <- B[, in_order, drop = FALSE]
    pi_in <- pi_in[in_order]
    block_sizes_in <- block_sizes_in[in_order]
  }

  pi_out <- pi_out / sum(pi_out)
  pi_in <- pi_in / sum(pi_in)

  # sample block memberships

  z_out <- sample_block_memberships(n, k_out, pi_out, block_sizes_out)
  z_out <- factor(z_out, levels = 1:k_out, labels = paste0(1:k_out))

  if (shared_blocks) {
    z_in <- z_out
  } else {
    z_in <- sample_block_memberships(n, k_in, pi_in, block_sizes_in)
    z_in <- factor(z_in, levels = 1:k_in, labels = paste0(1:k_in))
  }

  if (shared_blocks) {
    # outgoing and incoming parameters are tied together node by node, so
    # node i must correspond to row i of both X and Y. keep everything
    # aligned in node order, and permute all of it jointly when sorting.
    # (`shared_theta = TRUE` implies `shared_blocks = TRUE`)

    if (force_identifiability) {
      theta_out <- l1_normalize_within(theta_out, z_out)
      theta_in <- l1_normalize_within(theta_in, z_in)
    }

    if (sort_nodes) {
      node_order <- order(z_out, z_in, -theta_out, -theta_in)
      z_out <- z_out[node_order]
      z_in <- z_in[node_order]
      theta_out <- theta_out[node_order]
      theta_in <- theta_in[node_order]
    }

    if (k_out > 1) {
      X <- sparse.model.matrix(~ z_out + 0)
    } else {
      X <- Matrix(1, nrow = n, ncol = 1)
    }

    if (k_in > 1) {
      Y <- sparse.model.matrix(~ z_in + 0)
    } else {
      Y <- Matrix(1, nrow = n, ncol = 1)
    }

    X <- set_membership_weights(X, as.numeric(theta_out))
    Y <- set_membership_weights(Y, as.numeric(theta_in))

    dcsbm <- new_directed_dcsbm(
      X = X,
      S = B,
      Y = Y,
      theta_out = theta_out,
      theta_in = theta_in,
      z_out = z_out,
      z_in = z_in,
      pi_out = pi_out,
      pi_in = pi_in,
      sorted = sort_nodes,
      shared_blocks = shared_blocks,
      shared_theta = shared_theta,
      poisson_edges = poisson_edges,
      allow_self_loops = allow_self_loops,
      ...
    )

    return(validate_directed_dcsbm(dcsbm))
  }

  # co-factor model: outgoing and incoming parameters are unrelated, and
  # rows and columns are sorted independently

  if (sort_nodes) {
    z_out <- sort(z_out)
    z_in <- sort(z_in)
  }

  if (k_out > 1) {
    X <- sparse.model.matrix(~ z_out + 0)
  } else {
    X <- Matrix(1, nrow = n, ncol = 1)
  }

  if (k_in > 1) {
    Y <- sparse.model.matrix(~ z_in + 0)
  } else {
    Y <- Matrix(1, nrow = n, ncol = 1)
  }

  if (force_identifiability) {
    theta_out <- l1_normalize_within(theta_out, z_out)
    theta_in <- l1_normalize_within(theta_in, z_in)
  }

  X <- set_membership_weights(X, theta_out)
  Y <- set_membership_weights(Y, theta_in)

  if (sort_nodes) {
    # X, Y indexing must match z indexing, be careful about this
    X <- sort_by_all_columns(X)
    Y <- sort_by_all_columns(Y)

    if (k_out > 1) {
      theta_out <- X@x
    } else {
      theta_out <- sort(theta_out, decreasing = TRUE)
    }

    if (k_in > 1) {
      theta_in <- Y@x
    } else {
      theta_in <- sort(theta_in, decreasing = TRUE)
    }
  }

  dcsbm <- new_directed_dcsbm(
    X = X,
    S = B,
    Y = Y,
    theta_out = theta_out,
    theta_in = theta_in,
    z_out = z_out,
    z_in = z_in,
    pi_out = pi_out,
    pi_in = pi_in,
    sorted = sort_nodes,
    shared_blocks = shared_blocks,
    shared_theta = shared_theta,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops,
    ...
  )

  validate_directed_dcsbm(dcsbm)
}

#' @method print directed_dcsbm
#' @export
print.directed_dcsbm <- function(x, ...) {
  cat(glue("Directed Degree-Corrected Stochastic Blockmodel\n", .trim = FALSE))
  cat(glue(
    "-----------------------------------------------\n\n",
    .trim = FALSE
  ))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Outgoing Blocks (k_out): {x$k1}\n", .trim = FALSE))
  cat(glue("Incoming Blocks (k_in): {x$k2}\n", .trim = FALSE))

  shared_blocks <- if (isTRUE(x$shared_blocks)) {
    "TRUE (z_out and z_in are identical)"
  } else {
    "FALSE (z_out and z_in are independent)"
  }

  shared_theta <- if (isTRUE(x$shared_theta)) {
    "TRUE (theta_out and theta_in are identical)"
  } else {
    "FALSE (theta_out and theta_in are separate)"
  }

  cat(glue("Shared blocks: {shared_blocks}\n", .trim = FALSE))
  cat(glue("Shared degree heterogeneity: {shared_theta}\n\n", .trim = FALSE))

  cat("Traditional DCSBM parameterization:\n\n")
  cat("Block memberships (z_out):", dim_and_class(x$z_out), "\n")
  cat("Block memberships (z_in):", dim_and_class(x$z_in), "\n")
  cat("Degree heterogeneity (theta_out):", dim_and_class(x$theta_out), "\n")
  cat("Degree heterogeneity (theta_in):", dim_and_class(x$theta_in), "\n")
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
