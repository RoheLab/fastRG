library(igraph)

test_that("undirected factor model", {
  n <- 10000
  k <- 5

  X <- matrix(rpois(n = n * k, 1), nrow = n)
  S <- matrix(runif(n = k * k, 0, .1), nrow = k)

  ufm <- undirected_factor_model(X, S, expected_degree = 10)

  expect_equal(expected_degree(ufm), 10)

  mean(rowSums(sample_sparse(ufm)))

  expected_degree(ufm)

  graph <- sample_igraph(ufm)
  sample_mean_degree <- mean(degree(graph))
  sample_mean_degree

  expect_lt(9, sample_mean_degree)
  expect_lt(sample_mean_degree, 11)

  expect_silent(eigs_sym(ufm))
})

test_that("directed factor model", {

  n <- 10000
  d <- 1000

  k1 <- 5
  k2 <- 3

  X <- matrix(rpois(n = n * k1, 1), nrow = n)
  Y <- matrix(rpois(n = d * k2, 1), nrow = d)
  S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1)

  dfm <- directed_factor_model(X = X, S = S, Y = Y, expected_in_degree = 10)
  expect_equal(expected_in_degree(dfm), 10)

  dfm2 <- directed_factor_model(X = X, S = S, Y = Y, expected_out_degree = 10)
  expect_equal(expected_out_degree(dfm2), 10)

  dfm3 <- directed_factor_model(X = X, S = S, Y = Y, expected_density = 0.1)
  expect_equal(expected_density(dfm3), 0.1)

  expect_silent(svds(dfm))
})
