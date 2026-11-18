#' @importFrom RSpectra svds
#' @export
RSpectra::svds

#' @importFrom RSpectra eigs_sym
#' @export
RSpectra::eigs_sym

# below this size, RSpectra's Lanczos solvers are unreliable on the exactly
# low-rank operators that factor models produce, see
# https://github.com/yixuan/RSpectra/issues/20
small_model_threshold <- 25

rspectra_issue_url <- "https://github.com/yixuan/RSpectra/issues/20"

# RSpectra defaults to ncv = min(dim, max(2 * k + 1, 20)), which for small
# models means the Krylov subspace is (nearly) the whole space. on the exactly
# rank-k operators of factor models, this makes RSpectra error or silently
# return wrong values much more often. a smaller default is far more
# reliable. respect any ncv the user supplies
default_ncv <- function(opts, k, dim) {
  if (is.null(opts$ncv)) {
    ncv <- min(dim - 1, max(2 * k + 1, 10))

    # RSpectra requires k < ncv <= dim. when that isn't possible, leave the
    # choice to RSpectra
    if (ncv > k && ncv <= dim) {
      opts$ncv <- ncv
    }
  }

  opts
}

warn_if_small_model <- function(dim, call = rlang::caller_env()) {
  if (min(dim) <= small_model_threshold) {
    rlang::warn(
      c(
        paste0(
          "RSpectra can return inaccurate eigenvalues or singular values for ",
          "small factor models (here min(n, d) = ",
          min(dim),
          ")."
        ),
        i = paste0("See <", rspectra_issue_url, "> for details."),
        i = paste0(
          "For small models, compare with `eigen()` or `svd()` of ",
          "`as.matrix(expectation(model))`."
        )
      ),
      class = "fastRG_warning_small_model_spectra",
      call = call
    )
  }

  invisible(NULL)
}

check_no_sigma <- function(sigma, call = rlang::caller_env()) {
  if (!is.null(sigma)) {
    rlang::abort(
      c(
        "`sigma` (shift-and-invert mode) is not supported for factor models.",
        i = paste0(
          "Use `eigs_sym(as.matrix(expectation(model)), sigma = ...)` instead."
        )
      ),
      call = call
    )
  }

  invisible(NULL)
}

#' Compute the eigendecomposition of the expected adjacency matrix of an undirected factor model
#'
#' @param A An [undirected_factor_model()].
#' @param k Desired rank of decomposition.
#' @inherit RSpectra::eigs_sym params details
#' @param sigma Not supported for factor models, which are represented
#'   implicitly. Must be `NULL`. To use shift-and-invert mode, call
#'   `eigs_sym()` on `as.matrix(expectation(A))` instead.
#' @param opts Control parameters passed to [RSpectra::eigs_sym()]. See
#'   the "Accuracy for small models" section for the default `ncv`.
#' @param ... Passed to [RSpectra::eigs_sym()].
#'
#' @section Accuracy for small models:
#'
#' These methods compute spectra with RSpectra's Lanczos solvers, applied
#' implicitly to the expected adjacency matrix, which has rank `k`. RSpectra
#' is unreliable on such exactly low-rank matrices when they are small: it
#' can error with "TridiagEigen: eigen decomposition failed", or silently
#' return wrong eigenvalues or singular values (see
#' <https://github.com/yixuan/RSpectra/issues/20>).
#'
#' To reduce how often this happens, these methods default to a smaller
#' Krylov subspace than RSpectra does, `opts$ncv = min(dim - 1, max(2 * k + 1, 10))`,
#' unless you specify `opts$ncv` yourself. They also warn whenever
#' `min(n, d)` is 25 or less. For such small models, it is cheap and
#' reliable to compute the decomposition directly with [eigen()] or [svd()]
#' of `as.matrix(expectation(A))`.
#'
#' @method eigs_sym undirected_factor_model
#' @export
eigs_sym.undirected_factor_model <- function(
  A,
  k = A$k,
  which = "LM",
  sigma = NULL,
  opts = list(),
  ...
) {
  rlang::check_installed("RSpectra")
  check_no_sigma(sigma)
  warn_if_small_model(A$n)

  Ax <- function(x, args) as.numeric(args$X %*% (args$SXt %*% x))

  eigs_sym(
    Ax,
    k,
    which = which,
    opts = default_ncv(opts, k, A$n),
    ...,
    n = A$n,
    args = list(X = A$X, SXt = tcrossprod(A$S, A$X))
  )
}

#' Compute the singular value decomposition of the expected adjacency matrix of an undirected factor model
#'
#' @param A An [undirected_factor_model()].
#' @param k Desired rank of decomposition.
#' @inherit RSpectra::svds params details
#' @param opts Control parameters passed to [RSpectra::svds()]. See
#'   the "Accuracy for small models" section for the default `ncv`.
#' @param ... Passed to [RSpectra::svds()].
#' @inheritSection eigs_sym.undirected_factor_model Accuracy for small models
#'
#' @method svds undirected_factor_model
#' @export
svds.undirected_factor_model <- function(
  A,
  k = A$k,
  nu = k,
  nv = k,
  opts = list(),
  ...
) {
  rlang::check_installed("RSpectra")
  warn_if_small_model(A$n)

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
    opts = default_ncv(opts, k, A$n),
    ...,
    Atrans = Atx,
    dim = c(A$n, A$n),
    args = list(X = A$X, S = A$S)
  )
}

#' Compute the singular value decomposition of the expected adjacency matrix of a directed factor model
#'
#' @param A A [directed_factor_model()].
#' @param k Desired rank of decomposition.
#' @inherit RSpectra::svds params details
#' @param opts Control parameters passed to [RSpectra::svds()]. See
#'   the "Accuracy for small models" section for the default `ncv`.
#' @param ... Passed to [RSpectra::svds()].
#' @inheritSection eigs_sym.undirected_factor_model Accuracy for small models
#'
#' @method svds directed_factor_model
#' @export
svds.directed_factor_model <- function(
  A,
  k = min(A$k1, A$k2),
  nu = k,
  nv = k,
  opts = list(),
  ...
) {
  rlang::check_installed("RSpectra")
  warn_if_small_model(c(A$n, A$d))

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
    opts = default_ncv(opts, k, min(A$n, A$d)),
    ...,
    Atrans = Atx,
    dim = c(A$n, A$d),
    args = list(X = A$X, S = A$S, Y = A$Y)
  )
}
