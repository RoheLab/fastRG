#' Sample a directed stochastic blockmodel
#'
#' @param n Number of nodes in graph.
#'
#' @param pi_in Block sampling proportions for incoming blocks. Must be
#'   positive, but do not need to sum to one, as they will be
#'   normalized internally. Note that `length(pi_in)` implicitly
#'   specifies the number of incoming blocks. The number of
#'   incoming and outgoing blocks need not be the same.
#'
#' @param pi_out Block sampling proportions for outgoing blocks. Must be
#'   positive, but do not need to sum to one, as they will be
#'   normalized internally. Note that `length(pi_out)` implicitly
#'   specifies the number of outgoing blocks. The number of
#'   incoming and outgoing blocks need not be the same.
#'
#' @param B A `length(pi_in)` by `length(pi_out)` matrix of block connection
#'   probabilities. `B[i, j]` contains the probability that a node in
#'   community `i` links to a node in community `j`. Does not need to
#'   be symmetric.
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
#' @family directed graphs
#'
#' @details TODO: describe the generative model.
#'
#' @examples
#'
#' set.seed(27)
#'
#' n <- 1000
#' pi_in <- c(0.3, 0.7)
#' pi_out <- c(0.2, 0.6, 0.2)
#'
#' B <- rbind(
#'   c(0.7, 0.1, 0.3),
#'   c(0.1, 0.6, 0.2)
#' )
#'
#' A <- disbm(n, pi_in, pi_out, B)
#'
disbm <- function(n, pi_in, pi_out, B, avg_deg = NULL, poisson_edges = TRUE,
                  sort_nodes = FALSE, ...) {

  params <- disbm_params(
    n = n, pi_in = pi_in, pi_out = pi_out, B = B, avg_deg = avg_deg,
    poisson_edges = poisson_edges, sort_nodes = sort_nodes
  )

  # NOTE: avg_deg is null since we handled scaling B internally
  fastRG(params$X, params$S, params$Y, poisson_edges = poisson_edges,
         avg_deg = NULL, ...)
}

#' @export
disbm_params <- function(n, pi_in, pi_out, B, avg_deg = NULL,
                         poisson_edges = TRUE, sort_nodes = FALSE, ...) {

  if (any(pi_in < 0)) {
    stop("Elements of `pi_in` must be non-negative.", call. = FALSE)
  }

  if (any(pi_out < 0)) {
    stop("Elements of `pi_out` must be non-negative.", call. = FALSE)
  }

  pi_in <- pi_in / sum(pi_in)
  K_in <- length(pi_in)

  pi_out <- pi_out / sum(pi_out)
  K_out <- length(pi_out)

  if (K_in != nrow(B) || K_out != ncol(B)) {
    stop(
      "`B` must be a `length(pi_in)` x `length(pi_out)` matrix.",
      call. = FALSE
    )
  }

  # incoming block memberships
  z_in <- sample(K_in, n, replace = TRUE, prob = pi_in)
  z_out <- sample(K_out, n, replace = TRUE, prob = pi_out)

  z_in <- factor(z_in, levels = 1:K_in)
  z_out <- factor(z_out, levels = 1:K_out)

  # sort(z) orders the nodes so that all nodes in the first
  # block are together, nodes in the second block are all together, etc
  if (sort_nodes) {
    z_in <- sort(z_in)
    z_out <- sort(z_out)
  }

  X <- sparse.model.matrix(~ z_in + 0)
  Y <- sparse.model.matrix(~ z_out + 0)

  # now we handle avg_deg specially to make sure we don't get a matrix B
  # with probabilities outside of [0, 1]

  if (is.null(avg_deg)) {
    return(list(X = X, S = B, Y = Y))
  }

  # scale B just like in fastRG()
  B <- B * avg_deg / expected(X, B, Y)$degree

  if (!poisson_edges) {

    if (max(B) > 1)
      stop(
        "Expected edge values must be not exceed 1 for bernoulli graphs. ",
        "Either diminish `avg_deg` or set `poisson_edges = TRUE`.",
        call. = FALSE
      )

    # we're still sampling from a Poisson distribution, but B has been
    # specified as Bernoulli edges probabilities. convert these edges
    # probabilities such that we can feed them into a Poisson sampling
    # procedure

    B <- -log(1 - B)
  }

  list(X = X, S = B, Y = Y, z_in = z_in, z_out = z_out)
}
