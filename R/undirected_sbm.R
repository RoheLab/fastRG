
validate_undirected_sbm <- function(x) {

  values <- unclass(x)

  if (!inherits(x, "undirected_sbm")) {
    stop(
      "Undirected SBMs must inherit \"undirected_sbm\" class!",
      call. = FALSE
    )
  }

  if (!all(values$theta == 1)) {
    stop(
      "`theta` must equal one for all nodes in SBMs without degree ",
      "correction.",
      call. = FALSE
    )
  }

  if (!(values$edge_distribution %in% c("poisson", "bernoulli"))) {
    stop(
      "`edge_distribution` must be either \"poisson\" or \"bernoulli\".",
      call. = fALSE
    )
  }

  if (values$edge_distribution == "bernoulli" && max(values$S) > 1) {
    stop(
      "Elements of `B` must be not exceed 1 for bernoulli SBMs.",
      call. = FALSE
    )
  }

  x
}

#' Create an undirected stochastic blockmodel object
#'
#' To specify a stochastic blockmodel, you must specify
#' the number of nodes (via `n`), the mixing matrix (via `k` or `B`),
#' and the relative block probabilites (optional, via `pi`).
#' We provide sane defaults for most of these options to enable
#' rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n The number of nodes in the network. Must be
#'   a positive integer. This argument is required.
#'
#' @param edge_distribution Either `"poisson"` or `"bernoulli"`. The
#'   default is `"poisson"`, in which case the SBM can be a
#'   multigraph, i.e. multiple edges between the same two nodes
#'   are allowed. If `edge_distribution == "bernoulli"` only a
#'   single edge is allowed between any pair of nodes. See Section 2.3
#'   of Rohe et al (2017) for details.
#'
#' @inherit dcsbm params details
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#'
#' @references Rohe, Karl, Jun Tao, Xintian Han, and Norbert Binkiewicz. 2017.
#'    “A Note on Quickly Sampling a Sparse Matrix with Low Rank Expectation.”
#'    Journal of Machine Learning Research; 19(77):1-13, 2018.
#'    <http://www.jmlr.org/papers/v19/17-128.html>
#'
#' @return An `undirected_sbm` S3 object, which is a subclass of the
#'   [dcsbm()] object, with one additional field.
#'
#'   - `edge_distribution`: Either "poisson" or "bernoulli".
#'
#' @export
#' @seealso [fastRG()]
#' @family stochastic block models
#' @family undirected graphs
#'
#' @details
#'
#' A stochastic block is equivalent to a degree-corrected
#' stochastic blockmodel where the degree heterogeneity parameters
#' have all been set equal to 1.
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_sbm <- sbm(n = 1000, k = 5, expected_density = 0.01)
#' lazy_sbm
#'
#' # by default we get a multigraph (i.e. multiple edges are
#' # allowed between the same two nodes). using bernoulli edges
#' # will with an adjacency matrix with only zeroes and ones
#'
#' bernoulli_sbm <- sbm(
#'   n = 5000,
#'   k = 300,
#'   edge_distribution = "bernoulli",
#'   expected_degree = 80
#' )
#'
#' bernoulli_sbm
#'
#' edgelist <- sample_edgelist(bernoulli_sbm)
#' edgelist
#'
#' A <- sample_sparse(bernoulli_sbm)
#'
#' # only zeroes and ones!
#' unique(A)
#'
sbm <- function(
  n,
  k = NULL, B = NULL,
  ...,
  pi = rep(1 / k, k),
  edge_distribution = c("poisson", "bernoulli"),
  sort_nodes = TRUE) {

  edge_distribution <- rlang::arg_match(edge_distribution)

  sbm <- dcsbm(
    n = n,
    theta = rep(1, n),
    k = k,
    B = B,
    ...,
    pi = pi,
    sort_nodes = sort_nodes,
    subclass = "undirected_sbm"
  )

  sbm$edge_distribution <- edge_distribution

  if (edge_distribution == "bernoulli") {

    # we're still sampling from a Poisson distribution, but S has been
    # specified as Bernoulli edges probabilities. convert these edges
    # probabilities such that we can feed them into a Poisson sampling
    # procedure

    sbm$S <- -log(1 - sbm$S)
  }

  validate_undirected_sbm(sbm)
}

#' @method print undirected_sbm
#' @export
print.undirected_sbm <- function(x, ...) {

  cat(glue("Undirected Stochastic Blockmodel\n", .trim = FALSE))
  cat(glue("--------------------------------\n\n", .trim = FALSE))

  sorted <- if (x$sorted) "arranged by block" else "not arranged by block"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k}\n\n", .trim = FALSE))

  cat("Traditional SBM parameterization:\n\n")
  cat("Block memberships (z):", dim_and_class(x$z), "\n")
  cat("Block probabilities (pi):", dim_and_class(x$pi), "\n")
  cat(glue("Edge distribution: {x$edge_distribution}\n\n", .trim = FALSE))

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}


# dispatch hacks to respect type of edges---------------------------------------

#' @rdname sample_edgelist
#' @export
sample_edgelist.undirected_sbm <- function(
  factor_model,
  ...,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  poisson_edges <- factor_model$edge_distribution == "poisson"

  NextMethod()

  NextMethod("sample_edgelist", factor_model, ..., poisson_edges = FALSE)
}

#' @rdname sample_sparse
#' @export
sample_sparse.undirected_erdos_renyi <- function(
  factor_model,
  ...,
  poisson_edges = FALSE,
  allow_self_loops = TRUE) {

  NextMethod("sample_sparse", factor_model, ..., poisson_edges = FALSE)
}

#' @rdname sample_igraph
#' @export
sample_igraph.undirected_erdos_renyi <- function(
  factor_model,
  ...,
  poisson_edges = FALSE,
  allow_self_loops = TRUE) {

  NextMethod("sample_igraph", factor_model, ..., poisson_edges = FALSE)
}

#' @rdname sample_tidygraph
#' @export
sample_tidygraph.undirected_erdos_renyi <- function(
  factor_model,
  ...,
  poisson_edges = FALSE,
  allow_self_loops = TRUE) {

  NextMethod("sample_tidygraph", factor_model, ..., poisson_edges = FALSE)
}


#' @rdname sample_edgelist
#' @export
sample_edgelist.undirected_sbm <- function(
  factor_model,
  ...,
  allow_self_loops = TRUE) {

  poisson_edges <- factor_model$edge_distribution == "poisson"

  NextMethod()
}
