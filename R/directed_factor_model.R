new_directed_factor_model <- function(
  X, S, Y,
  poisson_edges,
  allow_self_loops,
  ...,
  subclass = character()) {

  ellipsis::check_dots_unnamed()

  n <- nrow(X)
  k1 <- ncol(X)
  d <- nrow(Y)
  k2 <- ncol(Y)

  model <- list(
    X = X,
    S = S,
    Y = Y,
    n = n,
    k1 = k1,
    d = d,
    k2 = k2,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )

  class(model) <- c(subclass, "directed_factor_model")
  model
}

validate_directed_factor_model <- function(x) {

  values <- unclass(x)

  if (any(values$X < 0) || any(values$S < 0) || any(values$Y < 0)) {
    stop(
      "`X`, `S` and `Y` can only contain non-negative elements.",
      call. = FALSE
    )
  }

  if (values$k1 != nrow(values$S)) {
    stop("`k1` must equal the number of rows in `S`", call. = FALSE)
  }

  if (values$k2 != ncol(values$S)) {
    stop("`k2` must equal the number of columns in `S`", call. = FALSE)
  }

  if (!is.logical(values$poisson_edges)) {
    stop("`poisson_edges` must be a logical(1) vector.", call. = FALSE)
  }

  if (!is.logical(values$allow_self_loops)) {
    stop("`allow_self_loops` must be a logical(1) vector.", call. = FALSE)
  }

  x
}

#' Create a directed factor model graph
#'
#' A directed factor model graph is a directed
#' generalized Poisson random dot product graph. The edges
#' in this graph are assumpted to be independent and Poisson
#' distributed. The graph is parameterized by its expected
#' adjacency matrix, with is `E[A] = X S Y'`. We do not recommend
#' that causal users use this function, see instead `directed_dcsbm()`
#' and related functions, which will formulate common variants
#' of the stochastic blockmodels as undirected factor models
#' *with lots of helpful input validation*.
#'
#' @param X A [matrix()] or [Matrix()] representing real-valued
#'   latent node positions encoding community structure of
#'   incoming edges. Entries must be positive.
#'
#' @param S A [matrix()] or [Matrix()] mixing matrix. Entries
#'   must be positive.
#'
#' @param Y A [matrix()] or [Matrix()] representing real-valued
#'   latent node positions encoding community structure of
#'   outgoing edges. Entries must be positive.
#'
#' @param ... Ignored. For internal developer use only.
#'
#' @param expected_in_degree If specified, the desired expected in degree
#'   of the graph. Specifying `expected_in_degree` simply rescales `S`
#'   to achieve this. Defaults to `NULL`. Specify only one of
#'   `expected_in_degree`, `expected_out_degree`, and `expected_density`.
#'
#' @param expected_out_degree If specified, the desired expected out degree
#'   of the graph. Specifying `expected_out_degree` simply rescales `S`
#'   to achieve this. Defaults to `NULL`. Specify only one of
#'   `expected_in_degree`, `expected_out_degree`, and `expected_density`.
#'
#' @param expected_density If specified, the desired expected density
#'   of the graph. Specifying `expected_density` simply rescales `S`
#'   to achieve this. Defaults to `NULL`. Specify only one of
#'   `expected_in_degree`, `expected_out_degree`, and `expected_density`.
#'
#' @param poisson_edges Logical indicating whether or not
#'   multiple edges are allowed to form between a pair of
#'   nodes. Defaults to `TRUE`. When `FALSE`, sampling proceeds
#'   as usual, and duplicate edges are removed afterwards. Further,
#'   when `FALSE`, we assume that `S` specifies a desired between-factor
#'   connection probability, and back-transform this `S` to the
#'   appropriate Poisson intensity parameter to approximate Bernoulli
#'   factor connection probabilities. See Section 2.3 of Rohe et al. (2017)
#'   for some additional details.
#'
#' @param allow_self_loops Logical indicating whether or not
#'   nodes should be allowed to form edges with themselves.
#'   Defaults to `TRUE`. When `FALSE`, sampling proceeds allowing
#'   self-loops, and these are then removed after the fact.
#'
#' @return A `directed_factor_model` S3 class based on a list
#'   with the following elements:
#'
#'   - `X`: The incoming latent positions as a [Matrix()] object.
#'
#'   - `S`: The mixing matrix as a [Matrix()] object.
#'
#'   - `Y`: The outgoing latent positions as a [Matrix()] object.
#'
#'   - `n`: The number of nodes with incoming edges in the network.
#'
#'   - `k1`: The dimension of the latent node position vectors
#'     encoding incoming latent communities (i.e. in `X`).
#'
#'   - `d`: The number of nodes with outgoing edges in the network.
#'     Does not need to match `n` -- rectangular adjacency matrices
#'     are supported.
#'
#'   - `k2`: The dimension of the latent node position vectors
#'     encoding outgoing latent communities (i.e. in `Y`).
#'
#'   - `poisson_edges`: Whether or not the graph is taken to be have
#'     Poisson or Bernoulli edges, as indicated by a logical vector
#'     of length 1.
#'
#'   - `allow_self_loops`: Whether or not self loops are allowed.
#'
#' @export
#'
#' @examples
#'
#' n <- 10000
#'
#' k1 <- 5
#' k2 <- 3
#'
#' d <- 5000
#'
#' X <- matrix(rpois(n = n * k1, 1), nrow = n)
#' S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
#' Y <- matrix(rexp(n = k2 * d, 1), nrow = d)
#'
#' fm <- directed_factor_model(X, S, Y)
#' fm
#'
#' fm2 <- directed_factor_model(X, S, Y, expected_in_degree = 50)
#' fm2
#'
directed_factor_model <- function(
  X, S, Y,
  ...,
  expected_in_degree = NULL,
  expected_out_degree = NULL,
  expected_density = NULL,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  X <- Matrix(X)
  S <- Matrix(S)
  Y <- Matrix(Y)

  degree_controls <- c(
    !is.null(expected_in_degree),
    !is.null(expected_out_degree),
    !is.null(expected_density)
  )

  if (sum(degree_controls) > 1) {
    stop(
      "Must specify only one of `expected_in_degree` ",
      "`expected_out_degree`, and `expected_density`.",
      call. = FALSE
    )
  }

  fm <- new_directed_factor_model(
    X, S, Y,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops,
    ...
  )

  if (!is.null(expected_in_degree)) {

    if (expected_in_degree <= 0) {
      stop(
        "`expected_in_degree` must be strictly greater than zero.",
        call. = FALSE
      )
    }

    S <- S * expected_in_degree / expected_in_degree(fm)
  }

  if (!is.null(expected_out_degree)) {

    if (expected_out_degree <= 0) {
      stop(
        "`expected_out_degree` must be strictly greater than zero.",
        call. = FALSE
      )
    }

    S <- S * expected_out_degree / expected_out_degree(fm)
  }

  if (!is.null(expected_density)) {

    if (expected_density <= 0 || 1 <= expected_density) {
      stop(
        "`expected_density` must be strictly between zero and one.",
        call. = FALSE
      )
    }

    S <- S * expected_density / expected_density(fm)
  }

  fm$S <- S

  if (!poisson_edges) {

    # when poisson_edges = FALSE, S is the desired Bernoulli edge probability.
    # we must
    # back-transform it to a Poisson parameterization of S. see section 2.3
    # of the paper and issue #20 for details.

    if (max(fm$S) > 1) {
      stop(
        "Elements of `S` (after symmetrizing and scaling to achieve expected ",
        "degree) must not exceed 1 for Bernoulli graphs.",
        call. = FALSE
      )
    }

    fm$S <- -log(1 - fm$S)
  }

  validate_directed_factor_model(fm)
}

dim_and_class <- function(x, ...) {

  if (is.matrix(x) || inherits(x, "Matrix"))
    paste0(nrow(x), " x ", ncol(x), " [", class(x)[1], "]")
  else
    paste0(length(x), " [", class(x)[1], "]")
}

#' @method print directed_factor_model
#' @export
print.directed_factor_model <- function(x, ...) {

  cat(glue("Directed Factor Model\n", .trim = FALSE))
  cat(glue("---------------------\n\n", .trim = FALSE))

  cat(glue("Incoming Nodes (n): {x$n}\n", .trim = FALSE))
  cat(glue("Incoming Rank (k1): {x$k1}\n", .trim = FALSE))
  cat(glue("Outgoing Rank (k2): {x$k2}\n", .trim = FALSE))
  cat(glue("Outgoing Nodes (d): {x$d}\n\n", .trim = FALSE))

  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n")
  cat("Y:", dim_and_class(x$Y), "\n\n")

  cat("Poisson edges:", as.character(x$poisson_edges), "\n")
  cat("Allow self loops:", as.character(x$allow_self_loops), "\n\n")

  cat(
    glue("Expected edges: {round(expected_edges(x))}"),
    glue("Expected density: {round(expected_density(x), 5)}"),
    glue("Expected in degree: {round(expected_in_degree(x), 1)}"),
    glue("Expected out degree: {round(expected_out_degree(x), 1)}"),
    sep = "\n"
  )
}

