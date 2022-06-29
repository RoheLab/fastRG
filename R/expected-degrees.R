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
  X <- factor_model$X
  S <- factor_model$S
  as.numeric(X %*% tcrossprod(S, X))
}

#' @export
expected_density.undirected_factor_model <- function(factor_model, ...) {
  expected_edges(factor_model) / as.numeric(choose(factor_model$n, 2))
}

#' @importFrom RSpectra eigs_sym
#' @export
RSpectra::eigs_sym

#' @method eigs_sym undirected_factor_model
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

#' @method svds undirected_factor_model
#' @export
svds.undirected_factor_model <- function(
    A,
    k = A$k,
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
    as.numeric(args$X %*% (tcrossprod(args$S, args$X) %*% x))
  }

  Atx <- function(x, args) {
    as.numeric(tcrossprod(args$X, args$S) %*% crossprod(args$X, x))
  }

  svds(
    A = Ax,
    k = k,
    nu = nu,
    nv = nv,
    opts = opts,
    ...,
    Atrans = Atx,
    dim = c(A$n, A$n),
    args = list(X = A$X, S = A$S)
  )
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

#' @importFrom RSpectra svds
#' @export
RSpectra::svds




#' @method svds directed_factor_model
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
    as.numeric(args$X %*% (tcrossprod(args$S, args$Y) %*% x))
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
