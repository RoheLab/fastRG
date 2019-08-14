#' Sample a degree corrected mixed membership stochastic blockmodel graph
#'
#' @param theta TODO
#' @param B TODO
#' @param alpha TODO
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
#'
#'   to remove DC, set theta=rep(1,n)
#'   xi  = theta_i * z_i, where z_i ~ dirichlet(alpha)
#'   lambda_ij = xi' B xj
#'   probability of i connecting to j:  1 - exp(-lambda_ij)
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
#' theta <- round(rlnorm(n, 1))  # degree parameter
#' alpha <- rep(1, k)            # prior parameter for nodes in communities
#'
#' A <- dc_mixed(theta, alpha, B, avgDeg = 50)
#'
dc_mixed <- function(theta, alpha, B, ...) {

  if (length(theta) == 1) {
    n <- theta
    theta <- rep(1, n)
  }

  if (length(theta) > 1) n <- length(theta)
  K <- length(alpha)

  if (K != nrow(B) || ncol(B) != nrow(B))
    stop("Both dimensions of B must match length of alpha", call. = FALSE)

  X <- t(igraph::sample_dirichlet(n, alpha))
  X <- X[order(X %*% (1:K)), ]
  Theta <- Diagonal(n, theta)
  X <- Theta %*% X

  fastRG(X, B, ...)
}
