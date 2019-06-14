#' Sample a degree corrected overlapping stochastic blockmodel graph
#'
#' @param theta TODO
#' @param pi TODO
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
#'   xi  = theta_i * z_i, where [z_i]_j ~ bernoulli(pi_j)
#'   lambda_ij = xi' B xj
#'   probability of i connecting to j:  1 - exp(-lambda_ij)
#'
#' @examples
#'
#' # TODO
#'
dcOverlapping <- function(theta, pi, B, ...) {

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
    X[, i] <- rbinom(n, 1, pi[i])

  X <- X[order(X %*% (1:K)), ]
  Theta <- Diagonal(n, theta)
  X <- Theta %*% X

  # TODO: deal with this irritating parameter stuff
  #   probably factor out into a different function
  # if (parametersOnly) return(list(X = X, S = B))

  # TODO: are there any arguments that don't make sense for this model
  # that should be ignored?
  fastRG(X, B, ...)
}


#' Sample a degree corrected mixed membership stochastic blockmodel graph
#'
#' @param theta TODO
#' @param pi TODO
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
#'
#'   to remove DC, set theta=rep(1,n)
#'   xi  = theta_i * z_i, where z_i ~ dirichlet(alpha)
#'   lambda_ij = xi' B xj
#'   probability of i connecting to j:  1 - exp(-lambda_ij)
#'
#' @examples
#'
#' # TODO
#'
dcMixed <- function(theta, alpha, B, ...) {

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


#' Sample a degree corrected stochastic blockmodel graph
#'
#' @param theta TODO
#' @param pi TODO
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
#'   samples from a degree corrected stochastic blockmodel
#'   pi a K vector, contains block sampling proportions
#'   theta an n vector, contains degree parameters
#'   B is K time K
#'
#'   Define lambda_ij = theta[i] * theta[j]*B_{U,V}
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
#' # TODO
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

#' Sample a stochastic blockmodel graph
#'
#' @param n TODO
#' @param pi Block sampling proportions.
#' @param B "Contains probabilities". TODO
#'
#' @inheritParams fastRG
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
#' # TODO
#'
sbm <- function(n, pi, B, avgDeg = NULL, PoissonEdges = TRUE,...) {

  K <- length(pi)

  if (K != nrow(B) || ncol(B) != nrow(B))
    stop("Both dimensions of B must match length of `pi`.", call. = FALSE)

  z <- sample(K, n, replace = TRUE, prob = pi)

  # you might want to comment this next line out...
  # but it is here so that pictures are pretty before clustering:
  z <- sort(z)
  X <- model.matrix(~ factor(as.character(z), levels = as.character(1:K)) - 1)

  # now we handle avgDeg specially to make sure we don't get a matrix B
  # with probabilities outside of [0, 1]

  if (is.null(avgDeg))
    return(fastRG(X, B))

  # scale B just like in fastRG()
  if (!is.null(avgDeg)) {
    eDbar <- howManyEdges(X, B)["avDeg"]
    B <- B * avgDeg / eDbar
  }


  if (!PoissonEdges && max(B) >= 1) {

    # TODO: probabilities of what?
    warning(
      "This combination of B and avgDeg has led to probabilities exceeding 1.",
      "We suggest you either diminish avgDeg or enable poisson edges."
    )

    # TODO: what is the following doing?
    # this ensure that edge probabilites are bijand not 1-exp(-bij)
    B <- -log(1 - B)
  }

  # NOTE: avgDeg is null since we handled scaling B internally
  fastRG(X, B, PoissonEdges = PoissonEdges, avgDeg = NULL, ...)
}


#' Sample an Erdos-Renyi graph
#'
#' @param n TODO
#' @param p TODO
#' @param avgDeg TODO
#' @param directed Defaults to `FALSE` for Erdos-Renyi graphs. The default
#'   in the more general [fastRG()] is `TRUE`.
#'
#' @inheritParams fastRG
#' @inherit fastRG params return
#'
#' @return Never returns Poisson edges.
#'
#' @export
#'
#' @seealso [fastRG()]
#' @family bernoulli graphs
#'
#' @examples
#'
#' # TODO
#'
er <- function(n, p = NULL, avgDeg = NULL, directed = FALSE, ...) {

  X <- matrix(1, nrow = n, ncol = 1)

  if (is.null(p) && is.null(avgDeg))
    stop("Must specify either `avgDeg` or `p`.", call. = FALSE)

  if (is.null(p))
    p <- avgDeg / n

  poisson_p <- -log(1 - p)
  S <- matrix(poisson_p, nrow = 1, ncol = 1)

  fastRG(X, S, PoissonEdges = FALSE, ...)
}

#' Sample a Chung-Lu graph
#'
#'
#' @param theta TODO
#'
#' @inheritParams fastRG
#' @inherit fastRG params return
#'
#' @return Never returns Poisson edges.
#'
#' @export
#'
#' @seealso [fastRG()]
#' @family bernoulli graphs
#'
#' @examples
#'
#' # TODO
#'
cl <- function(theta, ...) {
  # samples Chung-Lu graph
  # defaults to undirected.  If prefer directed, then define directed = TRUE
  # another

  # TODO: assumes theta is a vector and coerces to matrix to make fastRG()
  # happy?
  n <- length(theta)
  X <- matrix(theta, nrow = n, ncol = 1)
  S <- matrix(1, nrow = 1, ncol = 1)

  fastRG(X, S, PoissonEdges = FALSE, ...)
}
