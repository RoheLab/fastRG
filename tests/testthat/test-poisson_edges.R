library(igraph)
library(magrittr)

test_that("undirected graphs poisson_edges = FALSE", {
  set.seed(6)

  library(dplyr)

  n <- 1000
  k <- 5

  X <- matrix(rpois(n = n * k, 1), nrow = n)
  S <- matrix(runif(n = k * k, 0, .1), nrow = k)

  ufm <- undirected_factor_model(
    X, S,
    expected_density = 0.1, poisson_edges = FALSE
  )

  edgelist <- sample_edgelist(ufm)

  max_element_A <- edgelist %>%
    count(to, from, sort = TRUE) %>%
    pull(n) %>%
    max()

  expect_equal(max_element_A, 1)

  A <- sample_sparse(ufm)
  expect_equal(max(A), 1)

  igraph <- sample_igraph(ufm)
  expect_equal(max(as_adjacency_matrix(igraph)), 1)

  ### sampling graphs as tidygraph graphs ---------------

  tbl_graph <- sample_tidygraph(ufm)
  expect_equal(max(as_adjacency_matrix(tbl_graph)), 1)
})

test_that("directed graphs poisson_edges = FALSE", {
  set.seed(7)

  library(dplyr)

  n2 <- 1000

  k1 <- 5
  k2 <- 3

  d <- 500

  X <- matrix(rpois(n = n2 * k1, 1), nrow = n2)
  S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
  Y <- matrix(rexp(n = k2 * d, 1), nrow = d)

  fm <- directed_factor_model(X, S, Y, expected_density = 0.01, poisson_edges = FALSE)

  edgelist <- sample_edgelist(fm)

  max_element_A <- edgelist %>%
    count(to, from) %>%
    pull(n) %>%
    max()

  expect_equal(max_element_A, 1)

  A <- sample_sparse(fm)
  expect_equal(max(A), 1)

  igraph <- sample_igraph(fm)
  expect_equal(max(as_adjacency_matrix(igraph)), 1)

  ### sampling graphs as tidygraph graphs ---------------

  tbl_graph <- sample_tidygraph(fm)
  expect_equal(max(as_adjacency_matrix(tbl_graph)), 1)
})
