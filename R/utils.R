#' Calculate the expected number of edges in Poisson RDPG graph
#'
#' @param X An `n` by `k_1` matrix.
#' @param S A `k_1` by `k_2` matrix.
#' @param Y A `d` by `k_2` matrix. Defaults to `X`.
#'
#' @return A named list with three measurements from a Poisson RDPG(X, S, Y)
#'     multi-graph:
#'
#'   - `count`: Expected number of edges.
#'
#'   - `degree`: Expected degree of across all nodes TODO: for directed graphs
#'     is this the in-degree or the out-degree?
#'
#'   - `density`: Expected edge density. Equals the expected number of edges
#'     divided by the total number of edges.
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
expected <- function(X, S, Y = X) {

  if (any(X < 0) || any(S < 0) || any(Y < 0)) {
    stop("`X`, `S`, `Y` can only contain non-negative elements.", call. = FALSE)
  }

  Cx <- diag(colSums(X), nrow = ncol(X), ncol = ncol(X))
  Cy <- diag(colSums(Y), nrow = ncol(Y), ncol = ncol(Y))

  # convert integers to numerics to fend off integer overflow
  nx <- as.numeric(nrow(X))
  ny <- as.numeric(nrow(Y))

  edges <- sum(Cx %*% S %*% Cy)
  degree <- edges / nx

  # TODO: set `density = NULL` and issue a warning on overflows

  list(
    edges = edges,
    degree = degree,
    density = edges / (nx * ny)
  )
}

# TODO: ask Karl if there is a better way to do this
#' @rdname expected
#' @export
expected_degrees <- function(X, S, Y = X) {

  if (any(X < 0) || any(S < 0) || any(Y < 0)) {
    stop("`X`, `S`, `Y` can only contain non-negative elements.", call. = FALSE)
  }

  SYt <- tcrossprod(S, Y)
  in_degree <- as.numeric(X %*% rowSums(SYt))   # rowSums(A)
  out_degree <- as.numeric(colSums(X) %*% SYt)  # colSums(A)

  list(in_degree = in_degree, out_degree = out_degree)
}

# after the object orientation transformation, make this an
# S3 generic for RSpectra::svds()

#' @rdname expected
#' @export
expected_svds <- function(X, S, Y = X, ..., k = NULL) {

  if (any(X < 0) || any(S < 0) || any(Y < 0)) {
    stop("`X`, `S`, `Y` can only contain non-negative elements.", call. = FALSE)
  }

  n <- nrow(X)
  d <- nrow(Y)

  if (is.null(k))
    k <- ncol(X)

  # using closures rather than args

  SYt <- tcrossprod(S, Y)

  Ax <- function(x, args) X %*% (SYt %*% x)
  Atx <- function(x, args) (x %*% X) %*% SYt

  svds(Ax, k, Atrans = Atx, dim = c(n, d))
}

#' @rdname expected
#' @export
expected_eigs_sym <- function(X, S, ..., k = NULL) {

  if (any(X < 0) || any(S < 0)) {
    stop("`X`, `S`, `Y` can only contain non-negative elements.", call. = FALSE)
  }

  n <- nrow(X)

  if (is.null(k))
    k <- ncol(S)

  SXt <- tcrossprod(S, X)

  Ax <- function(x, args) as.numeric(args$X %*% (args$SXt %*% x))

  eigs_sym(Ax, k, n = n, args = list(X = X, SXt = SXt))
}

