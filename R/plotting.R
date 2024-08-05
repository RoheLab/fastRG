#' Calculate the expected adjacency matrix
#'
#' @param model A [directed_factor_model()] or an [undirected_factor_model()]
#'   object.
#' @param ... Unused.
#'
#' @return The expected value of the adjacency matrix, conditional on the
#'   latent factors `X` and `Y` (if the model is directed).
#' @export
#'
expectation <- function(model, ...) {
  UseMethod("expectation")
}

#' @rdname expectation
#' @export
expectation.undirected_factor_model <- function(model, ...) {
  model$X %*% tcrossprod(model$S, model$X)
}

#' @rdname expectation
#' @export
expectation.directed_factor_model <- function(model, ...) {
  model$X %*% tcrossprod(model$S, model$Y)
}

#' Plot (expected) adjacency matrices
#'
#' @inheritParams expectation
#' @param A A [matrix()], [Matrix::Matrix()] or [Matrix::sparseMatrix()] object.
#'
#' @return A [ggplot2::ggplot2()] plot.
#' @export
#'
#' @examples
#'
#' set.seed(27)
#'
#' model <- dcsbm(n = 10, k = 2, expected_density = 0.2)
#'
#' plot_expectation(model)
#'
#' A <- sample_sparse(model)
#'
#' plot_sparse_matrix(A)
#'
plot_expectation <- function(model) {
  EA <- as.matrix(expectation(model))
  if (is.null(rownames(EA))) {
    rownames(EA) <- 1:nrow(EA)
  }
  if (is.null(colnames(EA))) {
    colnames(EA) <- 1:ncol(EA)
  }
  plot_dense_matrix(EA)
}

#' @rdname plot_expectation
#' @import ggplot2
#' @export
plot_dense_matrix <- function(A, ...) {
  long <- dplyr::as_tibble(A, rownames = "row")
  long <- tidyr::gather(long, col, value, -row)
  long <- dplyr::mutate_all(long, as.numeric)

  ggplot(long, aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    theme_minimal() +
    labs(
      fill = "Expected edges"
    ) +
    theme_void()
}

#' @rdname plot_expectation
#' @export
plot_sparse_matrix <- function(A) {

  stopifnot(inherits(A, "sparseMatrix"))

  A <- methods::as(A, "CsparseMatrix")
  A <- methods::as(A, "generalMatrix")

  ggplot(summary(A), aes(x = i, y = j, fill = as.factor(x))) +
    geom_tile() +
    scale_fill_grey() +
    scale_y_reverse() +
    expand_limits(x = nrow(A), y = ncol(A)) +
    theme_void() +
    labs(
      fill = "Edges"
    )
}
