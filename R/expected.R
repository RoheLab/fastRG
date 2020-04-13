#' Calculate the expected number of edges in Poisson RDPG graph
#'
#' underestimate when edge distribution is bernoulli
#'
#' @param factor_model TODO
#' @param ... Ignored.
#'
#' @details Note that the running time of [fastRG()] is proportional to
#'   `count`, the expected number of edges in the graph.
#'
#' @export
#'
#' @examples
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
#' expected(X, S, Y)
#' expected_degrees(X, S, Y)
#' expected_svds(X, S, Y)
#'
expected_edges <- function(factor_model, ...) {
  ellipsis::check_dots_empty()
  UseMethod("expected_edges")
}

expected_degree <- function(factor_model, ...) {
  UseMethod("expected_degree")
}

expected_degrees <- function(factor_model, ...) {
  UseMethod("expected_degrees")
}

expected_in_degree <- function(factor_model, ...) {
  UseMethod("expected_in_degree")
}

expected_out_degree <- function(factor_model, ...) {
  UseMethod("expected_out_degree")
}

expected_density <- function(factor_model, ...) {
  UseMethod("expected_density")
}

#' @export
expected_edges.undirected_factor_model <- function(factor_model, ...) {

  X <- factor_model$X
  S <- factor_model$S

  Cx <- Diagonal(n = ncol(X), x = colSums(X))
  sum(Cx %*% S %*% Cx)
}

#' @export
expected_degree.undirected_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(factor_model$n)
}

#' @export
expected_degrees.undirected_factor_model <- function(factor_model, ...) {
  X <- factor_model$X
  S <- factor_model$S
  as.numeric(X %*% tcrossprod(S, X))
}

#' @export
expected_density.undirected_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(factor_model$n)^2
}

#' @importFrom RSpectra eigs_sym
#' @export
RSpectra::eigs_sym

#' @export
eigs_sym.undirected_factor_model <- function(
  A, k = A$k,
  which = "LM", sigma = NULL,
  opts = list(),
  ...) {

  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop(
      "Must install `RSpectra` for this functionality.",
      call. = FALSE
    )
  }

  Ax <- function(x, args) as.numeric(args$X %*% (args$SXt %*% x))

  eigs_sym(Ax, k, n = A$n, args = list(X = A$X, SXt = tcrossprod(A$S, A$X)))
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

# TODO: sanity check if dividing by n and d are right in the following

#' @export
expected_in_degree.directed_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(factor_model$n)
}

#' @export
expected_out_degree.directed_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(factor_model$d)
}

#' @export
expected_density.directed_factor_model <- function(factor_model, ...) {

  n <- factor_model$n
  d <- factor_model$d

  expected_edges(factor_model) / (as.numeric(n) * as.numeric(d))
}

#' @importFrom RSpectra svds
#' @export
RSpectra::svds

#' @export
svds.directed_factor_model <- function(
  A,
  k = min(A$k1, A$k2),
  nu = k,
  nv = k,
  opts = list(),
  ...) {

  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop(
      "Must install `RSpectra` for this functionality.",
      call. = FALSE
    )
  }

  Ax <- function(x, args) {
    as.numeric(args$X %*% (tcrossprod(S, Y) %*% x))
  }

  Atx <- function(x, args) {
    as.numeric(tcrossprod(args$Y, args$S) %*% crossprod(args$X, x))
  }

  svds(
    A = Ax,
    k = k,
    nu = nu,
    nv = nv,
    opts = opts,
    ...,
    Atrans = Atx,
    dim = c(A$n, A$d),
    args = list(X = A$X, S = A$S, Y = A$Y)
  )
}
