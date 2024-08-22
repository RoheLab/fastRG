#' Calculate the expected edges in Poisson RDPG graph
#'
#' These calculations are conditional on the latent factors
#' `X` and `Y`.
#'
#' @inherit sample_edgelist references params
#'
#' @return Expected edge counts, or graph densities.
#'
#' @details Note that the runtime of the `fastRG` algorithm is proportional to
#'   the expected number of edges in the graph. Expected edge count will be
#'   an underestimate of expected number of edges for Bernoulli
#'   graphs. See the Rohe et al for details.
#'
#' @export
#'
#' @examples
#'
#' ##### an undirected blockmodel example
#'
#' n <- 1000
#' pop <- n / 2
#' a <- .1
#' b <- .05
#'
#' B <- matrix(c(a, b, b, a), nrow = 2)
#'
#' b_model <- fastRG::sbm(n = n, k = 2, B = B, poisson_edges = FALSE)
#'
#' b_model
#'
#' A <- sample_sparse(b_model)
#'
#' # compare
#' mean(rowSums(triu(A)))
#'
#' pop * a + pop * b # analytical average degree
#'
#' ##### more generic examples
#'
#' n <- 10000
#' k <- 5
#'
#' X <- matrix(rpois(n = n * k, 1), nrow = n)
#' S <- matrix(runif(n = k * k, 0, .1), nrow = k)
#'
#' ufm <- undirected_factor_model(X, S)
#'
#' expected_edges(ufm)
#' expected_degree(ufm)
#' eigs_sym(ufm)
#'
#' n <- 10000
#' d <- 1000
#'
#' k1 <- 5
#' k2 <- 3
#'
#' X <- matrix(rpois(n = n * k1, 1), nrow = n)
#' Y <- matrix(rpois(n = d * k2, 1), nrow = d)
#' S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1)
#'
#' dfm <- directed_factor_model(X = X, S = S, Y = Y)
#'
#' expected_edges(dfm)
#' expected_in_degree(dfm)
#' expected_out_degree(dfm)
#'
#' svds(dfm)
#'
expected_edges <- function(factor_model, ...) {
  ellipsis::check_dots_empty()
  UseMethod("expected_edges")
}

#' @export
expected_edges.directed_factor_model <- function(factor_model, ...) {
  X <- factor_model$X
  S <- factor_model$S
  Y <- factor_model$Y

  Cx <- Diagonal(n = ncol(X), x = colSums(X))
  Cy <- Diagonal(n = ncol(Y), x = colSums(Y))
  sum(Cx %*% S %*% Cy)
}

#' @export
expected_edges.undirected_factor_model <- function(factor_model, ...) {
  X <- factor_model$X
  S <- factor_model$S

  Cx <- Diagonal(n = ncol(X), x = colSums(X))
  sum(Cx %*% S %*% Cx)
}

#' @rdname expected_edges
#' @export
expected_degree <- function(factor_model, ...) {
  UseMethod("expected_degree")
}


#' @rdname expected_edges
#' @export
expected_in_degree <- function(factor_model, ...) {
  UseMethod("expected_in_degree")
}

#' @rdname expected_edges
#' @export
expected_out_degree <- function(factor_model, ...) {
  UseMethod("expected_out_degree")
}

#' @rdname expected_edges
#' @export
expected_density <- function(factor_model, ...) {
  UseMethod("expected_density")
}

#' @export
expected_degree.undirected_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(factor_model$n)
}


#' @rdname expected_edges
#' @export
expected_degrees <- function(factor_model, ...) {
  UseMethod("expected_degrees")
}

#' @export
expected_degrees.undirected_factor_model <- function(factor_model, ...) {
  # rowSums of E[A|X, S] = XSX' are XSX'1 for 1 a column vector of ones
  # want to avoid memory cost of instantiating all of E[A|X, S], which is
  # typically large and dense
  X <- factor_model$X
  S <- factor_model$S
  ones <- matrix(1, nrow = nrow(X))
  as.numeric(X %*% (tcrossprod(S, X) %*% ones))
}

#' @export
expected_density.undirected_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(choose(factor_model$n, 2))
}

#' @export
expected_in_degree.directed_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(factor_model$d)
}

#' @export
expected_out_degree.directed_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(factor_model$n)
}

#' @export
expected_density.directed_factor_model <- function(factor_model, ...) {
  n <- factor_model$n
  d <- factor_model$d

  expected_edges(factor_model) / (as.numeric(n) * as.numeric(d))
}
