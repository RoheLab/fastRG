#' Sample a stochastic blockmodel graph
#'
#' @param n Number of nodes in graph.
#' @param pi Block sampling proportions. Must be positive, but do not need to
#'   sum to one, as they will be normalized internally. Note that `length(pi)`
#'   implicitly specifies the number of blocks.
#' @param B "Contains probabilities". TODO
#' @param avgDeg TODO
#' @param PoissonEdges TODO
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block. This incurs the cost of a fast
#'   integer radix sort, but even this may be prohibitive when sampling
#'   large graphs. Defaults to `TRUE`.
#'
#' @inheritDotParams fastRG
#' @inherit fastRG params return
#'
#' @export
#' @seealso [fastRG()]
#' @family stochastic block models
#'
#' @details
#'
#'   TODO
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
#'
sbm <- function(n, pi, B, avgDeg = NULL, PoissonEdges = TRUE,
                sort_nodes = TRUE, ...) {

  params <- sbm_params(n = n, pi = pi, B = B, avgDeg = avgDeg,
                       PoissonEdges = PoissonEdges, sort_nodes = sort_nodes)

  # NOTE: avgDeg is null since we handled scaling B internally
  fastRG(params$X, params$S, PoissonEdges = PoissonEdges, avgDeg = NULL, ...)
}

#' @rdname sbm
#' @export
sbm_params <- function(n, pi, B, avgDeg = NULL, PoissonEdges = TRUE,
                sort_nodes = TRUE, ...) {

  if (any(pi < 0))
    stop("Elements of `pi` must be non-negative.", call. = FALSE)

  pi <- pi / sum(pi)
  K <- length(pi)

  if (K != nrow(B) || ncol(B) != nrow(B))
    stop("`B` must be a kxk matrix, where k is `length(pi)`.", call. = FALSE)

  # block memberships
  z <- sample(K, n, replace = TRUE, prob = pi)

  # sort(z) orders the nodes so that all nodes in the first
  # block are together, nodes in the second block are all together, etc
  if (sort_nodes)
    z <- sort(z)

  # X is a dummy matrix for the block memberships
  X <- sparse.model.matrix(~ as.factor(z) - 1)

  # now we handle avgDeg specially to make sure we don't get a matrix B
  # with probabilities outside of [0, 1]

  if (is.null(avgDeg))
    return(list(X = X, S = B, Y = X))

  # scale B just like in fastRG()

  eDbar <- howManyEdges(X, B)[["avDeg"]]
  B <- B * avgDeg / eDbar

  if (!PoissonEdges && max(B) >= 1) {

    # TODO: probabilities of what?
    warning(
      "This combination of B and avgDeg has led to probabilities exceeding 1.",
      "We suggest you either diminish avgDeg or enable poisson edges."
    )

    # we're still sampling from a Poisson distribution, but the B has been
    # specified as Bernoulli edges probabilities. convert these edges
    # probabilities such that we can feed them into a Poisson sampling
    # procedure
    B <- -log(1 - B)
  }

  list(X = X, S = B, Y = X)
}

