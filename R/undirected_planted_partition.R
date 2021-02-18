validate_undirected_planted_partition <- function(x) {

  values <- unclass(x)

  if (!inherits(x, "undirected_planted_partition")) {
    stop(
      "Undirected planted partition models must inherit \"undirected_planted_partition\" class!",
      call. = FALSE
    )
  }

  if (x$k > x$n) {
    stop("`k` must be less than or equal to `n`.", call. = FALSE)
  }

  if (!is.null(x$within_block) && length(x$within_block) != 1) {
    stop("`within_block` must be a scalar.", call. = FALSE)
  }

  if (!is.null(x$between_block) && length(x$between_block) != 1) {
    stop("`between_block` must be a scalar.", call. = FALSE)
  }

  if (!is.null(x$within_block) && (x$within_block < 0 || x$within_block > 1)) {
    stop("`within_block` must be between zero and one.", call. = FALSE)
  }

  if (!is.null(x$between_block) && (x$between_block < 0 || x$between_block > 1)) {
    stop("`between_block` must be between zero and one.", call. = FALSE)
  }

  if (!is.null(x$a) && x$a < 0) {
    stop("`a` must be greater than zero.", call. = FALSE)
  }

  if (!is.null(x$b) && x$b < 0) {
    stop("`b` must be greater than zero.", call. = FALSE)
  }


  # diagonal B must be constant

  # off-diagonal of B must be constant

  x
}

#' Create an undirected planted partition object
#'
#' To specify a planted partition model, you must specify
#' the number of nodes (via `n`), the mixing matrix (optional, either via
#' `within_block/between_block` or `a/b`),
#' and the relative block probabilites (optional, via `pi`).
#' We provide sane defaults for most of these options to enable
#' rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param k Number of planted partitions, as a positive integer.
#'   This argument is required.
#'
#' @param within_block Probability of within block edges. Must be
#'   strictly between zero and one. Must specify either
#'   `within_block` and `between_block`, or `a` and `b` to determine
#'   edge probabilities.
#'
#' @param between_block Probability of between block edges. Must be
#'   strictly between zero and one. Must specify either
#'   `within_block` and `between_block`, or `a` and `b` to determine
#'   edge probabilities.
#'
#' @param a Integer such that `a/n` is the probability of edges
#'   within a block. Useful for sparse graphs. Must specify either
#'   `within_block` and `between_block`, or `a` and `b` to determine
#'   edge probabilities.
#'
#' @param b Integer such that `b/n` is the probability of edges
#'   between blocks. Useful for sparse graphs. Must specify either
#'   `within_block` and `between_block`, or `a` and `b` to determine
#'   edge probabilities.
#'
#' @inherit sbm params details references
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#'
#' @return An `undirected_planted_partition` S3 object, which is a subclass
#'   of the [sbm()] object, with additional fields:
#'
#'   - `within_block`: TODO
#'
#'   - `between_block`: TODO
#'
#' @export
#' @family stochastic block models
#' @family undirected graphs
#'
#' @details
#'
#' A planted partition model is stochastic blockmodel in which
#' the diagonal and the off-diagonal of the mixing matrix `B`
#' are both constant. This means that edge probabilities
#' depend only on whether two nodes belong to the same block,
#' or to different blocks, but the particular blocks themselves
#' don't have any impact apart from this.
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_pp <- planted_partition(
#'   n = 1000,
#'   k = 5,
#'   expected_density = 0.01,
#'   within_block = 0.1,
#'   between_block = 0.01
#' )
#'
#' lazy_pp
#'
planted_partition <- function(
  n,
  k,
  ...,
  within_block = NULL,
  between_block = NULL,
  a = NULL,
  b = NULL,
  pi = rep(1 / k, k),
  edge_distribution = c("poisson", "bernoulli"),
  sort_nodes = TRUE) {

  if (is.null(within_block)) {
    within_block <- a / n
  }

  if (is.null(between_block)) {
    between_block <- b / n
  }

  B <- matrix(between_block, nrow = k, ncol = k)
  diag(B) <- rep(within_block, k)

  pp <- sbm(
    n = n,
    k = k,
    B = B,
    ...,
    pi = pi,
    sort_nodes = sort_nodes,
    edge_distribution = edge_distribution,
    # subclass = "undirected_planted_partition"  # temporary hack
  )

  pp$within_block <- within_block
  pp$between_block <- between_block

  # validate_undirected_planted_partition(pp)
  pp
}

#' @method print undirected_planted_partition
#' @export
print.undirected_planted_partition <- function(x, ...) {

  cat(glue("Undirected Planted Partition Model\n", .trim = FALSE))
  cat(glue("----------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k}\n\n", .trim = FALSE))

  cat("Traditional planted partition parameterization:\n\n")
  cat("Within block edge probability:", dim_and_class(x$within_block), "\n")
  cat("Between block edge probability:", dim_and_class(x$between_block), "\n")
  cat("Block memberships (z):", dim_and_class(x$z), "\n")
  cat("Block probabilities (pi):", dim_and_class(x$pi), "\n")
  cat(glue("Edge distribution: {x$edge_distribution}\n\n", .trim = FALSE))

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
