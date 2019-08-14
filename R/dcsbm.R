#' Sample a degree corrected stochastic blockmodel graph
#'
#' @param theta Node degree distribution parameter.
#' @param pi vector of length k, must be normalized. to allow independent
#'   probabilities of membership in each block use the [dc_overlapping()]
#'   model
#' @param B TODO -- import from SBM
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
#'   samples from a degree corrected stochastic blockmodel
#'   pi a K vector, contains block sampling proportions
#'   theta an n vector, contains degree parameters
#'   B is K time K
#'
#'   Define \eqn{lambda_ij = theta[i] * theta[j]*B_{U,V}}
#'   where U,V are blockmemberships of i and j, sampled from multinomial(pi)
#'   i connects to j with probability 1- exp( -lambda_ij)
#'
#'   if lambda_ij is small, then 1- exp( -lambda_ij) approx lambda_ij
#'   if avgDeg is set, then B is scaled so that the expected average
#'   degree is avgDeg.
#'
#'   in fact, avgDeg is a slight upper bound that is good when
#'   the graph is sparse.
#'
#'   to make function easy to parameterize, this function does not allow
#'   for different degree distributions between blocks.
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
#' A <- dcsbm(theta, pi, B, avgDeg = 50)
#'
dcsbm <- function(theta, pi, B, ...) {

  n <- length(theta)
  K <- length(pi)
  B <- B[order(pi), ]
  B <- B[, order(pi)]
  pi <- sort(pi / sum(pi))

  if (K != nrow(B) || ncol(B) != nrow(B))
    stop("Both dimensions of B must match length of `pi`.", call. = FALSE)

  z <- sample(K, n, replace = TRUE, prob = pi)

  # you might want to comment this next line out...
  # but it is here so that pictures are pretty before clustering:
  z <- sort(z)
  X <- sparse.model.matrix(~ as.factor(z) - 1)
  ct <- c(0, cumsum(table(z)))

  for (i in 1:K) {
    theta[(ct[i] + 1):ct[i + 1]] <- -sort(-theta[(ct[i] + 1):ct[i + 1]])
  }

  X@x <- theta

  fastRG(X, B, ...)
}

#' @rdname dcsbm
#' @export
dcsbm_params <- function(theta, pi, B) {

  n <- length(theta)
  K <- length(pi)
  B <- B[order(pi), ]
  B <- B[, order(pi)]
  pi <- sort(pi / sum(pi))

  if (K != nrow(B) || ncol(B) != nrow(B))
    stop("Both dimensions of B must match length of `pi`.", call. = FALSE)

  z <- sample(K, n, replace = TRUE, prob = pi)

  # you might want to comment this next line out...
  # but it is here so that pictures are pretty before clustering:
  z <- sort(z)
  X <- sparse.model.matrix(~ as.factor(z) - 1)
  ct <- c(0, cumsum(table(z)))

  for (i in 1:K) {
    theta[(ct[i] + 1):ct[i + 1]] <- -sort(-theta[(ct[i] + 1):ct[i + 1]])
  }

  X@x <- theta

  list(X = X, B = B, Y = X)
}
