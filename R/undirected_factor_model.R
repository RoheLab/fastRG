new_undirected_factor_model <- function(
  X, S,
  ...,
  subclass = character()) {

  ellipsis::check_dots_unnamed()

  n <- nrow(X)
  k <- ncol(S)

  model <- list(
    X = X,
    S = S,
    n = n,
    k = k
  )

  class(model) <- c(subclass, "undirected_factor_model")
  model
}

validate_undirected_factor_model <- function(x) {

  values <- unclass(x)

  if (any(values$X < 0) || any(values$S < 0)) {
    stop(
      "`X` and `S` can only contain non-negative elements.",
      call. = FALSE
    )
  }

  if (values$n != nrow(values$X)) {
    stop("`n` must equal the number of rows in `X`", call. = FALSE)
  }

  if (values$k != ncol(values$S)) {
    stop("`k` must equal the number of columns in `S`", call. = FALSE)
  }

  if (ncol(values$X) != nrow(values$S)) {
    stop(
      "The number of columns of `X` must equal the number of rows of `S`.",
      call. = FALSE
    )
  }

  if (ncol(values$S) != ncol(values$X)) {
    stop(
      "The number of columns of `S` must equal the number of columns of `X`.",
      call. = FALSE
    )
  }

  if (ncol(values$S) != nrow(values$S)) {
    stop(
      "`S` must be square for undirected factor models.",
      call. = FALSE
    )
  }

  x
}

#' Create an undirected factor model graph
#'
#' An undirected factor model graph is an undirected
#' generalized Poisson random dot product graph. The edges
#' in this graph are assumed to be independent and Poisson
#' distributed. The graph is parameterized by its expected
#' adjacency matrix, which is `E[A|X] = X S X'`. We do not recommend
#' that casual users use this function, see instead [dcsbm()]
#' and related functions, which will formulate common variants
#' of the stochastic blockmodels as undirected factor models
#' *with lots of helpful input validation*.
#'
#' @param X A [matrix()] or [Matrix()] representing real-valued
#'   latent node positions. Entries must be positive.
#'
#' @param S A [matrix()] or [Matrix()] mixing matrix. `S` is
#'   symmetrized if it is not already, as this is the undirected
#'   case. Entries must be positive.
#'
#' @param ... Ignored. Must be empty.
#'
#' @param expected_degree If specified, the desired expected degree
#'   of the graph. Specifying `expected_degree` simply rescales `S`
#'   to achieve this. Defaults to `NULL`. Do not specify both
#'   `expected_degree` and `expected_density` at the same time.
#'
#' @param expected_density If specified, the desired expected density
#'   of the graph. Specifying `expected_density` simply rescales `S`
#'   to achieve this. Defaults to `NULL`. Do not specify both
#'   `expected_degree` and `expected_density` at the same time.
#'
#' @return An `undirected_factor_model` S3 class based on a list
#'   with the following elements:
#'
#'   - `X`: The latent positions as a [Matrix()] object.
#'
#'   - `S`: The mixing matrix as a [Matrix()] object.
#'
#'   - `n`: The number of nodes in the network.
#'
#'   - `k`: The rank of expectation matrix. Equivalently,
#'     the dimension of the latent node position vectors.
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
#' ufm
#'
#' sane_ufm <- undirected_factor_model(X, S, expected_degree = 50)
#' sane_ufm
#'
undirected_factor_model <- function(
  X, S,
  ...,
  expected_degree = NULL,
  expected_density = NULL) {

  X <- Matrix(X)
  S <- Matrix(S)

  if (!is.null(expected_degree) && !is.null(expected_density)) {
    stop(
      "Cannot specify both `expected_degree` and `expected_density`.",
      call. = FALSE
    )
  }

  ufm <- new_undirected_factor_model(X, S, ...)

  if (!is.null(expected_degree)) {

    if (expected_degree <= 0) {
      stop(
        "`expected_degree` must be strictly greater than zero.",
        call. = FALSE
      )
    }

    S <- S * expected_degree / expected_degree(ufm)
  }

  if (!is.null(expected_density)) {

    if (expected_density <= 0 || 1 <= expected_density) {
      stop(
        "`expected_density` must be strictly between zero and one.",
        call. = FALSE
      )
    }

    S <- S * expected_density / expected_density(ufm)
  }

  # symmetrize edge probabilities by setting S := (S + t(S)) / 2
  # then divide by 2 again this doubles edge probabilities. if S
  # is already symmetric, this leaves S unchanged

  S <- (S + t(S)) / 4

  ufm$S <- S

  validate_undirected_factor_model(ufm)
}

dim_and_class <- function(x, ...) {

  if (is.matrix(x) || inherits(x, "Matrix"))
    paste0(nrow(x), " x ", ncol(x), " [", class(x)[1], "]")
  else
    paste0(length(x), " [", class(x)[1], "]")
}

#' @method print undirected_factor_model
#' @export
print.undirected_factor_model <- function(x, ...) {

  cat(glue("Undirected Factor Model\n", .trim = FALSE))
  cat(glue("-----------------------\n\n", .trim = FALSE))

  cat(glue("Nodes (n): {x$n}\n", .trim = FALSE))
  cat(glue("Rank (k): {x$k}\n\n", .trim = FALSE))

  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}

