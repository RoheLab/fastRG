#' Sample a random dot product graph as a sparse Matrix
#'
#' @inherit sample_edgelist params details references examples description
#'
#' @return For undirected factor models, a sparse
#'   [Matrix::Matrix()] of class `dsCMatrix`. In particular,
#'   this means the `Matrix` object (1) has
#'   double data type, (2) is symmetric, and (3) is in
#'   column compressed storage format.
#'
#'   For directed factor models, a sparse
#'   [Matrix::Matrix()] of class `dgCMatrix`. This means
#'   the `Matrix` object (1) has double data type,
#'   (2) in *not* symmetric, and (3) is in column
#'   compressed storage format.
#'
#'   To reiterate: for undirected graphs, you will get
#'   a symmetric matrix. For directed graphs, you will
#'   get a general sparse matrix.
#'
#' @export
#' @family samplers
#'
sample_sparse <- function(
    factor_model,
    ...) {
  rlang::check_dots_unnamed()
  UseMethod("sample_sparse")
}

#' @rdname sample_sparse
#' @export
sample_sparse.undirected_factor_model <- function(
    factor_model,
    ...) {
  # to construct a symmetric sparseMatrix, we only pass in elements
  # of either the upper or lower diagonal (otherwise we'll get an error)
  # so we're going to sample as we want a directed graph, and then
  # flop all edges to the lower diagonal

  X <- factor_model$X
  S <- factor_model$S

  el <- sample_edgelist(factor_model, ...)
  n <- factor_model$n

  if (nrow(el) == 0) {
    return(sparseMatrix(1:n, 1:n, x = 0, dims = c(n, n)))
  }

  if (factor_model$poisson_edges) {
    # NOTE: x = 1 is correct to create a multigraph adjacency matrix
    # here. see ?Matrix::sparseMatrix for details, in particular the
    # documentation for arguments `i, j` and `x`
    #
    # in the poisson_edges = FALSE case, sample_edgelist will have
    # removed duplicate directed edges, but not duplicate undirected
    # edges, so we must still consider this case here
    A <- sparseMatrix(el$from, el$to, x = 1, dims = c(n, n), symmetric = TRUE)
  } else {
    A <- sparseMatrix(el$from, el$to, dims = c(n, n), symmetric = TRUE)
  }

  ### some comments on type-stability

  # there are three components of a Matrix object:
  #
  #   1. the data type (binary, integer, double)
  #   2. whether or not the Matrix is symmetric
  #   3. the storage format of the matrix
  #
  # the storage formats are: triplet, row compressed, column compressed
  #
  # sparseMatrix will give us a column-compressed (i.e. CSC) Matrix
  # by default, an a symmetric one, since we set symmetric = TRUE
  #
  # however the data type, and thus the corresponding Matrix class,
  # can vary. for bernoulli graphs (poisson_edges = FALSE) we
  # will get a binary data type, and for poisson graphs we get either
  # integer or double data (I forget which)
  #
  # further, Matrix provides tools to coerce between storage formats,
  # but casting between data types is only done implicitly. so
  # we are going to cast to matrix of doubles since that's
  # what has the widest support in extensions to the Matrix package

  A * 1.0 # LMAO
}

#' @rdname sample_sparse
#' @export
sample_sparse.directed_factor_model <- function(
    factor_model,
    ...) {
  edgelist <- sample_edgelist(factor_model, ...)

  n <- factor_model$n
  d <- factor_model$d

  if (nrow(edgelist) == 0) {
    return(sparseMatrix(1:n, 1:d, x = 0, dims = c(n, d)))
  }

  if (factor_model$poisson_edges) {
    # NOTE: x = 1 is correct to create a multigraph adjacency matrix
    # here. see ?Matrix::sparseMatrix for details, in particular the
    # documentation for arguments `i, j` and `x`

    A <- sparseMatrix(
      edgelist$from,
      edgelist$to,
      x = 1,
      dims = c(n, d),
      symmetric = FALSE
    )
  } else {
    A <- sparseMatrix(
      edgelist$from,
      edgelist$to,
      dims = c(n, d),
      symmetric = FALSE
    )
  }

  # see comments on type-stability of output in
  # sample_sparse.undirected_factor_model()
  #
  # here we will get double data type and CSC storage format
  # (which is as before), but the resulting matrix will no
  # longer be symmetric

  A * 1.0
}
