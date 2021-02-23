library(igraph)

test_that("undirected graphs are undirected", {

  set.seed(3)

  n <- 1000
  k <- 5

  X <- matrix(rpois(n = n * k, 1), nrow = n)
  S <- matrix(runif(n = k * k, 0, .1), nrow = k)

  ufm <- undirected_factor_model(
    X, S,
    expected_density = 0.1
  )

  edgelist <- sample_edgelist(ufm)
  expect_true(all(edgelist$from <= edgelist$to))

  A <- sample_sparse(ufm)
  expect_true(isSymmetric(A))

  igraph <- sample_igraph(ufm)
  expect_false(is_directed(igraph))

  ### sampling graphs as tidygraph graphs ---------------

  tbl_graph <- sample_tidygraph(ufm)
  expect_false(is_directed(tbl_graph))
})

test_that("directed graphs are directed", {

  set.seed(4)

  n2 <- 1000

  k1 <- 5
  k2 <- 3

  d <- 500

  X <- matrix(rpois(n = n2 * k1, 1), nrow = n2)
  S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
  Y <- matrix(rexp(n = k2 * d, 1), nrow = d)

  fm <- directed_factor_model(X, S, Y, expected_density = 0.01)

  edgelist <- sample_edgelist(fm)
  expect_false(all(edgelist$from <= edgelist$to))

  A <- sample_sparse(fm)
  expect_false(isSymmetric(A))

  igraph <- sample_igraph(fm)
  expect_false(is_directed(igraph))

  ### sampling graphs as tidygraph graphs ---------------

  tbl_graph <- sample_tidygraph(fm)
  expect_false(is_directed(tbl_graph))
})
