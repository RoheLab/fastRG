#' Sample a random dot product graph as an igraph graph
#'
#' There are two steps to using the `fastRG` package. First,
#' you must parameterize a random dot product graph by
#' specifying its expected adjacency matrix. Use functions such as
#' [dcsbm()], [sbm()], etc, to perform this specification.
#' Then, use [sample_igraph()] to generate a random graph,
#' represented as an [igraph::igraph()], with that expectation.
#'
#' @inherit sample_edgelist params details references examples
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

  graph_from_data_frame(edgelist, directed = FALSE)
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

  graph_from_data_frame(edgelist, directed = TRUE)
}
