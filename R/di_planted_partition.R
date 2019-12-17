#' Sample a directed, balanced planted partition model
#'
#' @param n Number of nodes in graph.
#'
#' @param k Number of incoming (and outgoing) planted partitions.
#'   There will be `k` incoming and `k` outgoing communities. In the future,
#'   we may allow different numbers of incoming and outgoing communities
#'   provided there is interest. Should be an integer.
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
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block. This incurs the the cost of a
#'   radix sort, which is linear in the number of nodes in the graph.
#'   Defaults to `FALSE`. Useful for plotting.
#'
#' @inheritDotParams fastRG
#'
#' @inherit fastRG params return
#'
#' @export
#' @seealso [fastRG()]
#' @family stochastic block models
#' @family bernoulli graphs
#' @family directed graphs
#'
#' @details TODO: describe the generative model.
#'
#' @examples
#'
#' set.seed(27)
#'
#' A <- di_planted_partition(
#'   n = 1000,
#'   k = 3,
#'   within_block = 0.8,
#'   between_block = 0.2
#' )
#'
#' B <- di_planted_partition(
#'   n = 1000,
#'   k = 3,
#'   a = 10,
#'   b = 4
#' )
#'
di_planted_partition <- function(n, k, within_block = NULL,
                                 between_block = NULL, a = NULL, b = NULL,
                                 sort_nodes = FALSE, ...) {

  params <- di_planted_partition_params(
    n = n, k = k, a = a, b = b, within_block = within_block,
    between_block = between_block, sort_nodes = sort_nodes
  )

  # NOTE: avg_deg is null since we handled scaling B internally
  fastRG(
    params$X, params$S, params$Y,
    poisson_edges = FALSE,
    directed = TRUE,
    ...
  )
}

#' @rdname planted_partition
#' @export
di_planted_partition_params <- function(n, k, within_block = NULL,
                                        between_block = NULL, a = NULL,
                                        b = NULL, sort_nodes = FALSE, ...) {

  if (k > n) {
    stop("`k` must be less than or equal to `n`.", call. = FALSE)
  }

  if (!is.null(within_block) && length(within_block) != 1) {
    stop("`within_block` must be a scalar.", call. = FALSE)
  }

  if (!is.null(between_block) && length(between_block) != 1) {
    stop("`between_block` must be a scalar.", call. = FALSE)
  }

  if (!is.null(within_block) && (within_block < 0 || within_block > 1)) {
    stop("`within_block` must be between zero and one.", call. = FALSE)
  }

  if (!is.null(between_block) && (between_block < 0 || between_block > 1)) {
    stop("`between_block` must be between zero and one.", call. = FALSE)
  }

  if (!is.null(a) && a < 0) {
    stop("`a` must be greater than zero.", call. = FALSE)
  }

  if (!is.null(b) && b < 0) {
    stop("`b` must be greater than zero.", call. = FALSE)
  }

  # TODO: more input validation on a and b

  if (is.null(within_block)) {
    within_block <- a / n
  }

  if (is.null(between_block)) {
    between_block <- b / n
  }

  # just a special case of the SBM, so leverage that infrastructure
  pi <- rep(1 / k, k)

  B <- matrix(between_block, nrow = k, ncol = k)
  diag(B) <- within_block

  disbm_params(n, pi, pi, B, poisson_edges = FALSE, sort_nodes = sort_nodes)
}
