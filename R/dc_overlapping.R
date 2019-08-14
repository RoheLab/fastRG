#' Sample a degree corrected overlapping stochastic blockmodel graph
#'
#' @param theta TODO
#' @param pi 1 by k vector. each element must be a probability.
#'   elements are considered independently, and do not need to sum to one
#' @param B TODO
#'
#' @inheritDotParams fastRG
#' @inherit fastRG return
#'
#' @export
#' @seealso [fastRG()]
#' @family stochastic block models
#'
#' @details
#'
#'   TODO: write the model out in detail here
#'
#'   to remove DC, set theta=rep(1,n)
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
#' B <- matrix(runif(k * k), nrow = k, ncol = k)  # mixing probabilities
#'
#' theta <- round(rlnorm(n, 2))  # degree parameter
#' pi <- c(1, 2, 4, 1, 1) / 5    # community membership
#'
#' A <- dc_overlapping(theta, pi, B, avgDeg = 50)
#'
dc_overlapping <- function(theta, pi, B, ...) {

  if (length(theta) == 1) {
    n <- theta
    theta <- rep(1, n)
  }

  if (length(theta) > 1)
    n <- length(theta)

  K <- length(pi)

  if (K != nrow(B) || ncol(B) != nrow(B))
    stop("Both dimensions of B must match length of `pi`.", call. = FALSE)

  X <- matrix(0, nrow = n, ncol = K)

  for (i in 1:K)
    X[, i] <- stats::rbinom(n, 1, pi[i])

  X <- X[order(X %*% (1:K)), ]
  Theta <- Diagonal(n, theta)
  X <- Theta %*% X

  fastRG(X, B, ...)
}
