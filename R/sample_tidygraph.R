#' Sample a random dot product graph as a tidygraph graph
#'
#' @inherit sample_edgelist params details references examples description
#'
#' @return A [tidygraph::tbl_graph()] object.
#'
#' @export
#' @family samplers
#'
sample_tidygraph <- function(
  factor_model,
  ...,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

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
  ...,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  edgelist <- sample_edgelist(
    factor_model,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )

  tidygraph::as_tbl_graph(edgelist, directed = FALSE)
}

#' @rdname sample_tidygraph
#' @export
sample_tidygraph.directed_factor_model <- function(
  factor_model,
  ...,
  poisson_edges = TRUE,
  allow_self_loops = TRUE) {

  edgelist <- sample_edgelist(
    factor_model,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )

  tidygraph::as_tbl_graph(edgelist, directed = TRUE)
}
