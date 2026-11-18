library(igraph)

test_that("undirected graphs are undirected", {
  set.seed(3)

  n <- 1000
  k <- 5

  X <- matrix(rpois(n = n * k, 1), nrow = n)
  S <- matrix(runif(n = k * k, 0, .1), nrow = k)

  ufm <- undirected_factor_model(
    X,
    S,
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

test_that("print.directed_factor_model() labels rows as outgoing nodes", {
  # rows of A send edges: n = 3 outgoing nodes, d = 5 incoming nodes
  m <- directed_factor_model(matrix(1, 3, 2), matrix(1, 2, 4), matrix(1, 5, 4))
  out <- paste(capture.output(print(m)), collapse = "\n")

  expect_match(out, "Outgoing Nodes (n): 3", fixed = TRUE)
  expect_match(out, "Outgoing Rank (k1): 2", fixed = TRUE)
  expect_match(out, "Incoming Rank (k2): 4", fixed = TRUE)
  expect_match(out, "Incoming Nodes (d): 5", fixed = TRUE)

  # out degree is edges per row, in degree is edges per column
  expect_equal(expected_out_degree(m), expected_edges(m) / 3)
  expect_equal(expected_in_degree(m), expected_edges(m) / 5)
})

test_that("directed_dcsbm() validator messages name the right side", {
  m <- directed_dcsbm(
    theta_out = rep(1, 10), theta_in = rep(1, 10), B = matrix(1, 2, 3)
  )

  bad <- m
  bad$pi_out <- c(0.5, 0.25, 0.25)
  expect_error(
    validate_directed_dcsbm(bad),
    "`pi_out` must match the outgoing rank",
    fixed = TRUE
  )

  bad <- m
  bad$pi_in <- c(0.5, 0.5)
  expect_error(
    validate_directed_dcsbm(bad),
    "`pi_in` must match the incoming rank",
    fixed = TRUE
  )
})
