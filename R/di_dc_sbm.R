#' Sample a directed degree corrected stochastic blockmodel graph
#'
#' @param theta_in A vector specifying the degree distribution parameters
#'   for node in-degree. Setting `theta = rep(1, n)` results in
#'   a stochastic blockmodel without in-degree correction.
#'
#' @param theta_out A vector specifying the degree distribution parameters
#'   for node out-degree. Setting `theta = rep(1, n)` results in
#'   a stochastic blockmodel without out-degree correction.
#'
#' @inherit disbm params
#' @inheritDotParams fastRG
#' @inherit fastRG return
#'
#' @export
#' @seealso [fastRG()]
#' @family stochastic block models
#' @family directed graphs
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
#' k_in <- 10
#' k_out <- 3
#'
#' B <- matrix(runif(k_in * k_out), nrow = k_in, ncol = k_out)
#'
#' theta_in <- round(rlnorm(n, 2))   # in degree heterogeneity
#' theta_out <- round(rlnorm(n, 1))  # out degree heterogeneity
#'
#' pi_in <- c(1, 1, 2, 3, 1, 3, 8, 2, 3, 1)
#' pi_out <- c(1, 2, 9)
#'
#' A <- di_dc_sbm(theta_in, theta_out, pi_in, pi_out, B, avg_deg = 50)
#'
#' params <- di_dc_sbm_params(theta_in, theta_out, pi_in, pi_out, B, avg_deg = 50)
#'
#' expected(params$X, params$S, params$Y)
#'
di_dc_sbm <- function(theta_in, theta_out, pi_in, pi_out, B, avg_deg = NULL,
                      poisson_edges = TRUE, sort_nodes = FALSE, ...) {

  params <- di_dc_sbm_params(
    pi_in = pi_in, pi_out = pi_out, theta_in = theta_in,
    theta_out = theta_out, B = B, avg_deg = avg_deg,
    poisson_edges = poisson_edges, sort_nodes = sort_nodes
  )

  # NOTE: avg_deg is null since we handled scaling B internally
  fastRG(params$X, params$S, params$Y, poisson_edges = poisson_edges,
         avg_deg = NULL, ...)
}

#' @export
di_dc_sbm_params <- function(theta_in, theta_out, pi_in, pi_out, B,
                             avg_deg = NULL, poisson_edges = TRUE,
                             sort_nodes = FALSE) {

  if (length(theta_in) != length(theta_out))
    stop(
      "Degree heterogeneity parameters must have same number of elements.",
      "Retry with `length(theta_in) == length(theta_out)`.",
      call. = FALSE
    )

  n <- length(theta_in)

  K_in <- length(pi_in)
  K_out <- length(pi_out)

  B <- B[order(pi_in), ]
  B <- B[, order(pi_out)]

  pi_in <- sort(pi_in / sum(pi_in))
  pi_out <- sort(pi_out / sum(pi_out))

  if (K_in != nrow(B) || K_out != ncol(B)) {
    stop(
      "`B` must be a `length(pi_in)` x `length(pi_out)` matrix.",
      call. = FALSE
    )
  }

  z_in <- sample(K_in, n, replace = TRUE, prob = pi_in)
  z_out <- sample(K_out, n, replace = TRUE, prob = pi_out)

  z_in <- factor(z_in, levels = 1:K_in)
  z_out <- factor(z_out, levels = 1:K_out)

  X <- sparse.model.matrix(~ z_in + 0)
  Y <- sparse.model.matrix(~ z_out + 0)

  if (sort_nodes) {
    z_in <- sort(z_in)
    z_out <- sort(z_out)

    # sort by degree within community as well as community

    ct_in <- c(0, cumsum(table(z_in)))

    for (i in 1:K_in) {
      theta_in[(ct_in[i] + 1):ct_in[i + 1]] <-
        -sort(-theta_in[(ct_in[i] + 1):ct_in[i + 1]])
    }

    ct_out <- c(0, cumsum(table(z_out)))

    for (i in 1:K_out) {
      theta_out[(ct_out[i] + 1):ct_out[i + 1]] <-
        -sort(-theta_out[(ct_out[i] + 1):ct_out[i + 1]])
    }

  }

  X@x <- theta_in
  Y@x <- theta_out

  if (is.null(avg_deg)) {
    return(list(X = X, S = B, Y = Y, z_in = z_in, z_out = z_out))
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
