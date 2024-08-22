#' @importFrom RSpectra svds
#' @export
RSpectra::svds

#' @importFrom RSpectra eigs_sym
#' @export
RSpectra::eigs_sym

#' Compute the eigendecomposition of the expected adjacency matrix of an undirected factor model
#'
#' @param A An [undirected_factor_model()].
#' @param k Desired rank of decomposition.
#' @inherit RSpectra::eigs_sym params details
#' @param ... Unused, included only for consistency with generic signature.
#'
#' @method eigs_sym undirected_factor_model
#' @export
eigs_sym.undirected_factor_model <- function(
    A, k = A$k,
    which = "LM", sigma = NULL,
    opts = list(),
    ...) {
  rlang::check_installed("RSpectra")

  Ax <- function(x, args) as.numeric(args$X %*% (args$SXt %*% x))

  eigs_sym(Ax, k, n = A$n, args = list(X = A$X, SXt = tcrossprod(A$S, A$X)))
}

#' Compute the singular value decomposition of the expected adjacency matrix of an undirected factor model
#'
#' @param A An [undirected_factor_model()].
#' @param k Desired rank of decomposition.
#' @inherit RSpectra::svds params details
#' @param ... Unused, included only for consistency with generic signature.
#'
#' @method svds undirected_factor_model
#' @export
svds.undirected_factor_model <- function(
    A,
    k = A$k,
    nu = k,
    nv = k,
    opts = list(),
    ...) {
  rlang::check_installed("RSpectra")

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

#' Compute the singular value decomposition of the expected adjacency matrix of a directed factor model
#'
#' @param A An [undirected_factor_model()].
#' @param k Desired rank of decomposition.
#' @inherit RSpectra::svds params details
#' @param ... Unused, included only for consistency with generic signature.
#'
#' @method svds directed_factor_model
#' @export
svds.directed_factor_model <- function(
    A,
    k = min(A$k1, A$k2),
    nu = k,
    nv = k,
    opts = list(),
    ...) {
  rlang::check_installed("RSpectra")

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
