#' Sample a balanced planted partition model
#'
#' @param n Number of nodes in graph.
#' @param k Number of planted partitions.
#' @param within_block Probability of within block edges.
#' @param between_block Probability of between block edges.
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
#'
#' @details TODO: describe the generative model.
#'
#' @examples
#'
#' set.seed(27)
#'
#' A <- planted_partition(
#'   n = 1000,
#'   k = 3,
#'   within_block = 0.8,
#'   between_block = 0.2
#' )
#'
planted_partition <- function(n, k, within_block, between_block,
                              sort_nodes = FALSE, ...) {

  params <- planted_partition_params(
    n = n, k = k, within_block = within_block,
    between_block = between_block, sort_nodes = sort_nodes
  )

  # NOTE: avg_deg is null since we handled scaling B internally
  fastRG(params$X, params$S, poisson_edges = FALSE, ...)
}

#' @rdname planted_partition
#' @export
planted_partition_params <- function(n, k, within_block, between_block,
                                     sort_nodes = FALSE, ...) {

  if (k > n) {
    stop("`k` must be less than or equal to `n`.", call. = FALSE)
  }

  if (length(within_block) != 1) {
    stop("`within_block` must be a scalar.", call. = FALSE)
  }

  if (length(between_block) != 1) {
    stop("`between_block` must be a scalar.", call. = FALSE)
  }

  if (within_block < 0 || within_block > 1) {
    stop("`within_block` must be between zero and one.", call. = FALSE)
  }

  if (between_block < 0 || between_block > 1) {
    stop("`between_block` must be between zero and one.", call. = FALSE)
  }

  # just a special case of the SBM, so leverage that infrastructure
  pi <- rep(1 / k, k)

  B <- matrix(between_block, nrow = k, ncol = k)
  diag(B) <- within_block

  sbm_params(n, pi, B, poisson_edges = FALSE, sort_nodes = sort_nodes)
}
