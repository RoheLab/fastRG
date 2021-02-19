#' Sample a random dot product graph as an igraph graph
#'
#' @inherit sample_edgelist params details references examples description
#'
#' @return An [igraph::igraph()] graph object.
#'
#' @export
#' @family samplers
#'
sample_igraph <- function(
  factor_model,
  ...,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  ellipsis::check_dots_unnamed()

  if (!(requireNamespace("igraph", quietly = TRUE))) {
    stop(
      "Must install `igraph` package to return graphs as `igraph` ",
      "objects",
      call. = FALSE
    )
  }

  UseMethod("sample_igraph")
}

#' @rdname sample_igraph
#' @export
sample_igraph.undirected_factor_model <- function(
  factor_model,
  ...,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  edgelist <- sample_edgelist(
    factor_model,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )

  igraph::graph_from_data_frame(edgelist, directed = FALSE)
}

#' @rdname sample_igraph
#' @export
sample_igraph.directed_factor_model <- function(
  factor_model,
  ...,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  edgelist <- sample_edgelist(
    factor_model,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )

  igraph::graph_from_data_frame(edgelist, directed = TRUE)
}
