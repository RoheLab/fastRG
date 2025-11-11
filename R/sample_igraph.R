#' Sample a random dot product graph as an igraph graph
#'
#' @inherit sample_edgelist params details references examples description
#'
#' @return An [igraph::igraph()] object that is possibly a
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
#'  Note that working with bipartite graphs in `igraph` is more
#'  complex than working with one-mode graphs.
#'
#' @export
#' @family samplers
#'
sample_igraph <- function(
  factor_model,
  ...
) {
  rlang::check_dots_unnamed()

  rlang::check_installed("igraph", "to return graphs as `igraph` objects.")

  UseMethod("sample_igraph")
}

#' @rdname sample_igraph
#' @export
sample_igraph.undirected_factor_model <- function(
  factor_model,
  ...
) {
  edgelist <- sample_edgelist(factor_model, ...)
  nodes <- tibble(name = 1:nrow(factor_model$X))
  igraph::graph_from_data_frame(edgelist, directed = FALSE, vertices = nodes)
}

#' @rdname sample_igraph
#' @export
sample_igraph.directed_factor_model <- function(
  factor_model,
  ...
) {
  if (factor_model$n == factor_model$d) {
    edgelist <- sample_edgelist(factor_model, ...)
    nodes <- tibble(name = 1:nrow(factor_model$X))
    one_mode <- igraph::graph_from_data_frame(
      edgelist,
      directed = TRUE,
      vertices = nodes
    )
    return(one_mode)
  }

  # to represent a rectangular adjacency matrix igraph we have to make
  # a bipartite graph rather than using the usual one-mode tools

  A <- sample_sparse(factor_model, ...)
  igraph::graph_from_biadjacency_matrix(A, directed = FALSE, multiple = TRUE)
}
