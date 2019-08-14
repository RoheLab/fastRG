#' Sample a stochastic blockmodel
#'
#' @param n Number of nodes in graph.
#'
#' @param pi Block sampling proportions. Must be positive, but do not need to
#'   sum to one, as they will be normalized internally. Note that `length(pi)`
#'   implicitly specifies the number of blocks.
#'
#' @param B A `length(pi)` by `length(pi)` matrix of block connection
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
#'
#' @details TODO: describe the generative model.
#'
#' @examples
#'
#' set.seed(27)
#'
#' n <- 1000
#' pi <- c(0.3, 0.7)
#'
#' B <- rbind(
#'   c(0.7, 0.1),
#'   c(0.1, 0.6)
#' )
#'
#' A <- sbm(n, pi, B)
sbm <- function(n, pi, B, avg_deg = NULL, poisson_edges = TRUE,
                sort_nodes = FALSE, ...) {
  params <- sbm_params(
    n = n, pi = pi, B = B, avg_deg = avg_deg,
    poisson_edges = poisson_edges, sort_nodes = sort_nodes
  )

  # NOTE: avg_deg is null since we handled scaling B internally
  fastRG(params$X, params$S, poisson_edges = poisson_edges, avg_deg = NULL, ...)
}

#' @rdname sbm
#' @export
sbm_params <- function(n, pi, B, avg_deg = NULL, poisson_edges = TRUE,
                       sort_nodes = FALSE, ...) {
  if (any(pi < 0)) {
    stop("Elements of `pi` must be non-negative.", call. = FALSE)
  }

  pi <- pi / sum(pi)
  K <- length(pi)

  if (K != nrow(B) || ncol(B) != nrow(B)) {
    stop("`B` must be a kxk matrix, where k is `length(pi)`.", call. = FALSE)
  }

  # block memberships
  z <- sample(K, n, replace = TRUE, prob = pi)

  # sort(z) orders the nodes so that all nodes in the first
  # block are together, nodes in the second block are all together, etc
  if (sort_nodes) {
    z <- sort(z)
  }

  # X is a dummy matrix for the block memberships
  X <- sparse.model.matrix(~ as.factor(z) - 1)

  # now we handle avg_deg specially to make sure we don't get a matrix B
  # with probabilities outside of [0, 1]

  if (is.null(avg_deg)) {
    return(list(X = X, S = B, Y = X))
  }

  # scale B just like in fastRG()
  B <- B * avg_deg / expected(X, B)$degree

  if (!poisson_edges && max(B) >= 1) {
    warning(
      "This combination of B and avg_deg has led to probabilities exceeding 1.",
      "We suggest you either diminish avg_deg or enable poisson edges.",
      call. = FALSE
    )

    # we're still sampling from a Poisson distribution, but the B has been
    # specified as Bernoulli edges probabilities. convert these edges
    # probabilities such that we can feed them into a Poisson sampling
    # procedure

    B <- -log(1 - B)
  }

  list(X = X, S = B, Y = X)
}
