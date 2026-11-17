#' Sample a random dot product graph as a tidygraph graph
#'
#' @inherit sample_edgelist params details references examples description
#'
#' @return A [tidygraph::tbl_graph()] object. For Poisson graphs
#'   (`poisson_edges = TRUE`) with `multigraph = TRUE` (the default), this
#'   is possibly a multigraph (that is, we take there to be multiple edges
#'   rather than weighted edges). With `multigraph = FALSE`, repeated edges
#'   are instead collapsed into a single edge, and the edge data has an
#'   integer `weight` column equal to the corresponding entry of the
#'   [sample_sparse()] adjacency matrix (all weights are `1` when
#'   `poisson_edges = FALSE`).
#'
#'   When `factor_model` is **undirected**:
#'
#'     - the graph is undirected and one-mode.
#'
#'   When `factor_model` is **directed** and **square**:
#'
#'     - the graph is directed and one-mode.
#'
#'   When `factor_model` is **directed** and **rectangular**:
#'
#'     - the graph is undirected and bipartite.
#'
#'  Note that working with bipartite graphs in `tidygraph` is more
#'  complex than working with one-mode graphs.
#'
#' @export
#' @family samplers
#'
sample_tidygraph <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  rlang::check_dots_unnamed()
  check_multigraph(multigraph)
  rlang::check_installed("igraph", "to return graphs as `tidygraph` objects.")

  UseMethod("sample_tidygraph")
}

#' @rdname sample_tidygraph
#' @export
sample_tidygraph.undirected_factor_model <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  ig <- sample_igraph(factor_model, ..., multigraph = multigraph)
  tidygraph::as_tbl_graph(ig, directed = FALSE)
}

#' @rdname sample_tidygraph
#' @export
sample_tidygraph.directed_factor_model <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  ig <- sample_igraph(factor_model, ..., multigraph = multigraph)
  tidygraph::as_tbl_graph(ig, directed = TRUE)
}
