#' Sample a degree corrected overlapping stochastic blockmodel graph
#'
#' @param p Community membership parameter. For each node, membership in
#'   each community is consider separately. Each node has probability
#'   `p[i]` of being in the `i^th` community. `length(p)` implicitly
#'   specifies the number of communities.
#'
#' @inherit dcsbm params
#' @inheritDotParams fastRG
#' @inherit fastRG return
#'
#' @export
#' @seealso [fastRG()]
#' @family stochastic block models
#'
#' @details TODO: write the model out in detail here
#'
#'   \deqn{
#'     xi  = \theta_i * z_i, where [z_i]_j ~ bernoulli(pi_j)
#'   }
#'
#'   \deqn{
#'     \lambda_{ij} = xi' B xj
#'   }
#'
#'   probability of \eqn{i} connecting to \eqn{j}:  \eqn{1 - exp(-\lambda_{ij})}
#'
#' @examples
#'
#' set.seed(27)
#'
#' n <- 1000
#' k <- 5
#'
#' B <- matrix(runif(k * k), nrow = k, ncol = k) # mixing probabilities
#'
#' theta <- round(rlnorm(n, 2))  # degree parameter
#' p <- c(1, 2, 4, 1, 1) / 5     # community membership
#'
#' A <- dc_overlapping(theta, p, B, avg_deg = 50)
#'
dc_overlapping <- function(theta, p, B, avg_deg = NULL,
                           poisson_edges = TRUE, sort_nodes = FALSE, ...) {

  params <- dc_overlapping_params(
    theta = theta, p = p, B = B,
    avg_deg = avg_deg,
    poisson_edges = poisson_edges,
    sort_nodes = sort_nodes
  )

  # NOTE: avg_deg is null since we handled scaling B internally
  fastRG(params$X, params$S, poisson_edges = poisson_edges, avg_deg = NULL, ...)
}

#' @rdname dc_overlapping
#' @export
dc_overlapping_params <- function(theta, p, B, avg_deg = NULL,
                                  poisson_edges = TRUE, sort_nodes = FALSE) {
  if (length(theta) == 1) {
    n <- theta
    theta <- rep(1, n)
  }

  if (length(theta) > 1) {
    n <- length(theta)
  }

  K <- length(p)

  if (K != nrow(B) || ncol(B) != nrow(B)) {
    stop("Both dimensions of B must match length of `pi`.", call. = FALSE)
  }

  X <- matrix(0, nrow = n, ncol = K)

  for (i in 1:K) {
    X[, i] <- stats::rbinom(n, 1, p[i])
  }

  if (sort_nodes)
    X <- X[order(X %*% (1:K)), ]

  Theta <- Diagonal(n, theta)
  X <- Theta %*% X

  list(X = X, S = B, Y = X)
}
