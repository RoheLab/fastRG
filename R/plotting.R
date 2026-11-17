#' Calculate the expected adjacency matrix
#'
#' @param model A [directed_factor_model()] or an [undirected_factor_model()]
#'   object.
#' @param ... Unused.
#'
#' @return The matrix \eqn{P = X S X^T}{P = X S X'} (undirected models) or
#'   \eqn{P = X S Y^T}{P = X S Y'} (directed models), as a
#'   [Matrix::Matrix()]. This is on the **Poisson scale**: it is the
#'   expected adjacency matrix, conditional on the latent factors, when
#'   edges are Poisson (`poisson_edges = TRUE`, the default), and not when
#'   edges are Bernoulli. See "Bernoulli graphs" below.
#'
#' @section Bernoulli graphs:
#'
#' When `poisson_edges = FALSE`, the model stores `S` on the Poisson scale:
#' the requested Bernoulli edge probabilities are transformed to Poisson
#' rates via \eqn{-\log(1 - p)}{-log(1 - p)} when the model is created.
#' Sampling then draws a Poisson number of edges for each pair of nodes and
#' keeps a single edge whenever that number is positive, so the probability
#' of an edge between nodes `i` and `j` is
#'
#' \deqn{
#'   \mathbb{E}[A_{ij} \mid X, S] = P(A_{ij} = 1 \mid X, S) = 1 - \exp(-P_{ij}),
#' }{
#'   E[A[i, j] | X, S] = P(A[i, j] = 1 | X, S) = 1 - exp(-P[i, j]),
#' }
#'
#' where \eqn{P}{P} is the matrix returned by `expectation()`. So for
#' Bernoulli graphs, compute the expected adjacency matrix as
#' `1 - exp(-as.matrix(expectation(model)))`. Since
#' \eqn{1 - e^{-x} \le x}{1 - exp(-x) <= x}, `expectation()` overestimates
#' edge probabilities for Bernoulli graphs, and the difference is small
#' when edge probabilities are small. For blockmodels, where each node
#' belongs to a single block with no degree heterogeneity, the result
#' recovers the requested probabilities exactly. If
#' `allow_self_loops = FALSE`, the diagonal of the expected adjacency matrix
#' is zero, whatever the edge distribution.
#'
#' The same Poisson scale applies to [expected_edges()],
#' [expected_degree()], [expected_density()] and the population spectra
#' computed by [eigs_sym()] and [svds()].
#'
#' @export
#'
#' @examples
#'
#' set.seed(27)
#'
#' B <- matrix(c(0.6, 0.1, 0.1, 0.4), nrow = 2)
#'
#' # Poisson edges: expectation() is E[A | X, S]
#' poisson_sbm <- sbm(n = 6, B = B, block_sizes = c(3, 3))
#' expectation(poisson_sbm)
#'
#' # Bernoulli edges: expectation() is on the Poisson scale
#' bernoulli_sbm <- sbm(n = 6, B = B, block_sizes = c(3, 3), poisson_edges = FALSE)
#' expectation(bernoulli_sbm)
#'
#' # the expected adjacency matrix, which recovers B for a blockmodel
#' 1 - exp(-as.matrix(expectation(bernoulli_sbm)))
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
  plot_dense_matrix(EA)
}

#' @rdname plot_expectation
#' @import ggplot2
#' @export
plot_dense_matrix <- function(A, ...) {
  # don't preserve rownames because rownames and colnames get coerced to
  # indices and this can lead to type errors
  rownames(A) <- 1:nrow(A)
  colnames(A) <- 1:ncol(A)

  long <- dplyr::as_tibble(A, rownames = "row")
  long <- tidyr::gather(long, col, value, -row)
  long <- dplyr::mutate_all(long, as.numeric)

  ggplot(long, aes(x = col, y = row, fill = value)) +
    geom_tile() +
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
