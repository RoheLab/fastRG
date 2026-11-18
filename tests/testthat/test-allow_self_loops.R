library(igraph)

test_that("undirected graphs allow_self_loops = FALSE", {
  set.seed(1)

  n <- 1000
  k <- 5

  X <- matrix(rpois(n = n * k, 1), nrow = n)
  S <- matrix(runif(n = k * k, 0, .1), nrow = k)

  ufm <- undirected_factor_model(
    X,
    S,
    expected_density = 0.1,
    allow_self_loops = FALSE
  )

  edgelist <- sample_edgelist(ufm, )
  expect_false(any(edgelist$from == edgelist$to))

  A <- sample_sparse(ufm)
  expect_false(any(diag(A) > 0))

  igraph <- sample_igraph(ufm)
  expect_false(any(diag(as_adjacency_matrix(igraph)) > 0))

  ### sampling graphs as tidygraph graphs ---------------

  tbl_graph <- sample_tidygraph(ufm)
  expect_false(any(diag(as_adjacency_matrix(tbl_graph)) > 0))
})

test_that("directed graphs allow_self_loops = FALSE", {
  set.seed(2)

  n2 <- 1000

  k1 <- 5
  k2 <- 3

  # square, so that row i and column i are the same node. rectangular
  # graphs have no self-loops, see the test below
  d <- 1000

  X <- matrix(rpois(n = n2 * k1, 1), nrow = n2)
  S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
  Y <- matrix(rexp(n = k2 * d, 1), nrow = d)

  fm <- directed_factor_model(
    X,
    S,
    Y,
    expected_density = 0.01,
    allow_self_loops = FALSE
  )

  edgelist <- sample_edgelist(fm)
  expect_false(any(edgelist$from == edgelist$to))

  A <- sample_sparse(fm)
  expect_false(any(diag(A) > 0))

  igraph <- sample_igraph(fm)
  expect_false(any(diag(as_adjacency_matrix(igraph)) > 0))

  ### sampling graphs as tidygraph graphs ---------------

  tbl_graph <- sample_tidygraph(fm)
  expect_false(any(diag(as_adjacency_matrix(tbl_graph)) > 0))
})

test_that("allow_self_loops = FALSE keeps all edges in rectangular models", {
  # rows and columns are different node sets, so from == to is a real edge
  set.seed(1)
  m <- directed_factor_model(
    matrix(1, 3, 1), matrix(2), matrix(1, 5, 1),
    allow_self_loops = FALSE
  )

  reps <- 500
  edges <- vapply(seq_len(reps), function(i) nrow(sample_edgelist(m)), numeric(1))
  diagonal <- vapply(seq_len(reps), function(i) {
    el <- sample_edgelist(m)
    sum(el$from == el$to)
  }, numeric(1))

  # total edges are Poisson(30); z-score of the mean
  mu <- expected_edges(m)
  expect_equal(mu, 30)
  expect_lt(abs(mean(edges) - mu) / sqrt(mu / reps), 5)

  # the 3 "diagonal" pairs each have rate 2
  expect_gt(sum(diagonal), 0)
  expect_lt(abs(mean(diagonal) - 6) / sqrt(6 / reps), 5)

  # square models still drop self-loops
  sq <- directed_factor_model(
    matrix(1, 3, 1), matrix(2), matrix(1, 3, 1),
    allow_self_loops = FALSE
  )
  el <- sample_edgelist(sq)
  expect_false(any(el$from == el$to))
})
