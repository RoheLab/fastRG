

#' @importFrom RSpectra svds
#' @export
RSpectra::svds


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
