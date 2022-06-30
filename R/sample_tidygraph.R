#' Sample a random dot product graph as a tidygraph graph
#'
#' @inherit sample_edgelist params details references examples description
#'
#' @return A [tidygraph::tbl_graph()] object that is possibly a
#'   multigraph (that is, we take there to be multiple edges
#'   rather than weighted edges).
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
  ...) {

  ellipsis::check_dots_unnamed()

  if (!(requireNamespace("tidygraph", quietly = TRUE))) {
    stop(
      "Must install `tidygraph` package to return graphs as `tidygraph` ",
      "objects",
      call. = FALSE
    )
  }

  UseMethod("sample_tidygraph")
}

#' @rdname sample_tidygraph
#' @export
sample_tidygraph.undirected_factor_model <- function(
  factor_model,
  ...) {

  nodes <- tibble(
    name = 1:nrow(factor_model$X)
  )

  edges <- sample_edgelist(factor_model, ...)
  tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
}

#' @rdname sample_tidygraph
#' @export
sample_tidygraph.directed_factor_model <- function(
  factor_model,
  ...) {

  ig <- sample_igraph(factor_model, ...)
  tidygraph::as_tbl_graph(ig, directed = TRUE)
}
