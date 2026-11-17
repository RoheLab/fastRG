#' Calculate the expected edges in Poisson RDPG graph
#'
#' These calculations are conditional on the latent factors
#' `X` and `Y`.
#'
#' @inherit sample_edgelist references params
#'
#' @return Expected edge counts, degrees, or graph densities.
#'
#' @details Write \eqn{P = \mathbb{E}[A \mid X, S] = X S X^T}{P = E[A | X, S] = X S X'}
#'   for an undirected factor model, and
#'   \eqn{P = X S Y^T}{P = X S Y'} for a directed one.
#'
#'   For **undirected** graphs, every edge between distinct nodes `i` and
#'   `j` appears twice in the symmetric adjacency matrix `A`, as both
#'   `A[i, j]` and `A[j, i]`, while a self-loop appears once, as `A[i, i]`.
#'   Accordingly:
#'
#'   - `expected_edges()` is the expected number of edges, counting each
#'     undirected edge once: \eqn{(\sum_{ij} P_{ij} + \mathrm{tr}(P)) / 2}{(sum(P) + trace(P)) / 2}.
#'     This includes expected self-loops, even when
#'     `allow_self_loops = FALSE`.
#'   - `expected_degree()` is the expected mean degree, i.e. the average
#'     row sum of \eqn{P}{P}: \eqn{\sum_{ij} P_{ij} / n}{sum(P) / n}.
#'     `expected_degrees()` gives the expected degree of each node.
#'   - `expected_density()` is the expected number of edges between
#'     distinct nodes, divided by the number of pairs of distinct nodes:
#'     \eqn{(\sum_{ij} P_{ij} - \mathrm{tr}(P)) / (n (n - 1))}{(sum(P) - trace(P)) / (n (n - 1))}.
#'     This is the average edge probability (or rate) between distinct
#'     nodes, and does not depend on self-loops.
#'
#'   For **directed** graphs, `expected_edges()` is
#'   \eqn{\sum_{ij} P_{ij}}{sum(P)}, `expected_out_degree()` and
#'   `expected_in_degree()` divide this by the number of rows and columns
#'   of `A`, and `expected_density()` divides it by the number of entries
#'   of `A`.
#'
#'   These are expectations under Poisson edges. For Bernoulli graphs
#'   (`poisson_edges = FALSE`), they are computed from the Poisson
#'   parameterization of `S`, and overestimate edge counts slightly when
#'   edge probabilities are large.
#'
#'   Note that the runtime of the `fastRG` algorithm is proportional to
#'   the expected number of edges in the graph. See Rohe et al. (2018) for
#'   details.
#'
#' @export
#'
#' @examples
#'
#' ##### an undirected blockmodel example
#'
#' n <- 100
#' pop <- n / 2
#' a <- .1
#' b <- .05
#'
#' B <- matrix(c(a, b, b, a), nrow = 2)
#'
#' b_model <- sbm(n = n, B = B, poisson_edges = FALSE)
#'
#' b_model
#'
#' A <- sample_sparse(b_model)
#'
#' # compare
#' mean(rowSums(A))
#'
#' pop * a + pop * b # analytical average degree
#'
#' ##### more generic examples
#'
#' n <- 1000
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
#' n <- 1000
#' d <- 100
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
  rlang::check_dots_empty()
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
  # each edge between distinct nodes is counted twice in sum(P), once as
  # P[i, j] and once as P[j, i], but self-loops are only counted once (#19)
  (expected_sum(factor_model) + expected_trace(factor_model)) / 2
}

# sum of all entries of P = E[A | X, S] = X S X', without forming P
expected_sum <- function(factor_model) {
  cx <- colSums(factor_model$X)
  as.numeric(crossprod(cx, factor_model$S %*% cx))
}

# trace of P = E[A | X, S] = X S X', i.e. the expected number of
# self-loops, without forming P
expected_trace <- function(factor_model) {
  X <- factor_model$X
  sum((X %*% factor_model$S) * X)
}

#' @rdname expected_edges
#' @export
expected_degrees <- function(factor_model, ...) {
  UseMethod("expected_degrees")
}

#' @export
expected_degrees.undirected_factor_model <- function(factor_model, ...) {
  # rowSums of E[A|X, S] = XSX' are XSX'1 for 1 a column vector of ones
  # want to avoid memory cost of instantiating all of E[A|X, S], which is
  # typically large and dense
  X <- factor_model$X
  S <- factor_model$S
  ones <- matrix(1, nrow = nrow(X))
  as.numeric(X %*% (tcrossprod(S, X) %*% ones))
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
  # mean row sum of P, i.e. the mean expected degree
  expected_sum(factor_model) / as.numeric(factor_model$n)
}

#' @export
expected_density.undirected_factor_model <- function(factor_model, ...) {
  # mean of P[i, j] over pairs of distinct nodes i != j
  n <- as.numeric(factor_model$n)
  off_diagonal <- expected_sum(factor_model) - expected_trace(factor_model)
  off_diagonal / (n * (n - 1))
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
