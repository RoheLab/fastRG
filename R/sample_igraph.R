#' Sample a random dot product graph as an igraph graph
#'
#' @inherit sample_edgelist params details references examples description
#'
#' @return An [igraph::igraph()] object. For Poisson graphs
#'   (`poisson_edges = TRUE`) with `multigraph = TRUE` (the default), this
#'   is possibly a multigraph (that is, we take there to be multiple edges
#'   rather than weighted edges). With `multigraph = FALSE`, repeated edges
#'   are instead collapsed into a single edge, and every edge has an integer
#'   `weight` attribute equal to the corresponding entry of the
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
#'  Note that working with bipartite graphs in `igraph` is more
#'  complex than working with one-mode graphs.
#'
#' @export
#' @family samplers
#'
sample_igraph <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  rlang::check_dots_unnamed()
  check_multigraph(multigraph)

  rlang::check_installed("igraph", "to return graphs as `igraph` objects.")

  UseMethod("sample_igraph")
}

#' @rdname sample_igraph
#' @export
sample_igraph.undirected_factor_model <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  edgelist <- sample_edgelist(factor_model, ..., multigraph = multigraph)
  nodes <- tibble(name = 1:nrow(factor_model$X))
  igraph::graph_from_data_frame(edgelist, directed = FALSE, vertices = nodes)
}

#' @rdname sample_igraph
#' @export
sample_igraph.directed_factor_model <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  if (factor_model$n == factor_model$d) {
    edgelist <- sample_edgelist(factor_model, ..., multigraph = multigraph)
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

  if (multigraph) {
    return(
      igraph::graph_from_biadjacency_matrix(A, directed = FALSE, multiple = TRUE)
    )
  }

  bipartite <- igraph::graph_from_biadjacency_matrix(
    A,
    directed = FALSE,
    weighted = TRUE
  )

  # entries of A are counts stored as doubles, keep weights integer to match
  # sample_edgelist(multigraph = FALSE)
  igraph::E(bipartite)$weight <- as.integer(igraph::E(bipartite)$weight)
  bipartite
}
