#' Sample a degree corrected stochastic blockmodel graph
#'
#' @param theta A vector specified the degree distribution parameters.
#'   The resulting graph will have `length(theta)` nodes. Setting
#'   `theta = rep(1, n)` recovers a stochastic blockmodel without
#'   degree correction.
#'
#' @inherit sbm params
#' @inheritDotParams fastRG
#' @inherit fastRG return
#'
#' @export
#' @seealso [fastRG()]
#' @family stochastic block models
#'
#' @details TODO: document the generative model
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
#' B <- matrix(runif(k * k), nrow = k, ncol = k) # mixing probabilities
#'
#' theta <- round(rlnorm(n, 2)) # degree parameter
#' pi <- c(1, 2, 4, 1, 1) / 5 # community membership
#'
#' A <- dcsbm(theta, pi, B, avg_deg = 50)
#'
#' params <- dcsbm_params(theta, pi, B, avg_deg = 50)
#'
#' expected(params$X, params$S)
#'
dcsbm <- function(theta, pi, B, avg_deg = NULL, poisson_edges = TRUE,
                  sort_nodes = FALSE, ...) {

  params <- dcsbm_params(
    theta = theta, pi = pi, B = B, avg_deg = avg_deg,
    poisson_edges = poisson_edges, sort_nodes = sort_nodes
  )

  # NOTE: avg_deg is null since we handled scaling B internally
  fastRG(params$X, params$S, poisson_edges = poisson_edges, avg_deg = NULL, ...)
}

#' @rdname dcsbm
#' @export
dcsbm_params <- function(theta, pi, B, avg_deg = NULL, poisson_edges = TRUE,
                         sort_nodes = FALSE) {

  # TODO: mimic the parameter rescaling from sbm()

  n <- length(theta)
  K <- length(pi)
  B <- B[order(pi), ]
  B <- B[, order(pi)]
  pi <- sort(pi / sum(pi))

  if (K != nrow(B) || ncol(B) != nrow(B)) {
    stop("Both dimensions of B must match length of `pi`.", call. = FALSE)
  }

  z <- sample(K, n, replace = TRUE, prob = pi)

  if (sort_nodes) {
    z <- sort(z)
  }

  X <- sparse.model.matrix(~ as.factor(z) - 1)

  # sort by degree within community as well as community

  if (sort_nodes) {
    ct <- c(0, cumsum(table(z)))

    for (i in 1:K) {
      theta[(ct[i] + 1):ct[i + 1]] <- -sort(-theta[(ct[i] + 1):ct[i + 1]])
    }

  }

  X@x <- theta

  if (is.null(avg_deg)) {
    return(list(X = X, S = B, Y = X))
  }

  # scale B just like in fastRG()
  B <- B * avg_deg / expected(X, B)$degree

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

  list(X = X, S = B, Y = X)
}
