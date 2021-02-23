test_that("undirected factor model", {

  set.seed(7)

  n <- 10000
  k <- 5

  # don't allow self edges at all in these calculations via
  # SBM model with zero on the diagonal of B

  B <- matrix(data = 0.5, nrow = k, ncol = k)
  diag(B) <- 0

  ufm <- sbm(n = n, k = k, B = B, expected_degree = 10)
  expect_equal(expected_degree(ufm), 10)
  expect_equal(expected_density(ufm), 0.001)


  el <- sample_edgelist(ufm)
  el_mean_degree <- 2 * nrow(el) / n
  expect_lt(9, el_mean_degree)
  expect_lt(el_mean_degree, 11)


  A <- sample_sparse(ufm)
  matrix_mean_degree <- mean(rowSums(A))

  expect_equal(rowSums(A), colSums(A))
  expect_lt(9, matrix_mean_degree)
  expect_lt(matrix_mean_degree, 11)


  graph <- sample_igraph(ufm)
  igraph_mean_degree <- mean(igraph::degree(graph))

  expect_lt(9, igraph_mean_degree)
  expect_lt(igraph_mean_degree, 11)


  tbl_graph <- sample_tidygraph(ufm)
  tbl_graph_mean_degree <- mean(igraph::degree(tbl_graph))

  expect_lt(9, tbl_graph_mean_degree)
  expect_lt(tbl_graph_mean_degree, 11)


  expect_silent(eigs_sym(ufm))
})

test_that("directed factor model", {

  set.seed(8)

  library(dplyr)

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


  ### edgelist tests -----------------------------------------------------------

  el <- sample_edgelist(dfm)

  el_mean_in_degree <- el %>%
    count(to) %>%
    pull(n) %>%
    mean()

  expect_lt(9, el_mean_in_degree)
  expect_lt(el_mean_in_degree, 11)



  el2 <- sample_edgelist(dfm2)

  el2_mean_out_degree <- el2 %>%
    count(from) %>%
    pull(n) %>%
    mean()

  expect_lt(95, el2_mean_out_degree)
  expect_lt(el2_mean_out_degree, 105)



  el3 <- sample_edgelist(dfm3)

  el3_density <- nrow(el3) / as.numeric(n * n)

  expect_lt(0.008, el3_density)
  expect_lt(el3_density, 0.012)


  ### sparse matrix tests ------------------------------------------------------

  A <- sample_sparse(dfm)
  matrix_mean_in_degree <- mean(colSums(A))

  expect_lt(9, matrix_mean_in_degree)
  expect_lt(matrix_mean_in_degree, 11)


  A2 <- sample_sparse(dfm2)
  matrix_mean_out_degree <- mean(rowSums(A2))

  expect_lt(95, matrix_mean_out_degree)
  expect_lt(matrix_mean_out_degree, 105)


  A3 <- sample_sparse(dfm3)
  A3_density <- mean(A3)

  expect_lt(0.008, A3_density)
  expect_lt(A3_density, 0.012)

  ### igraph tests --------------------------------------------------------------

  ig <- sample_igraph(dfm)
  ig_mean_in_degree <- mean(igraph::degree(ig, mode = "in"))

  expect_lt(9, ig_mean_in_degree)
  expect_lt(ig_mean_in_degree, 11)


  ig2 <- sample_igraph(dfm2)
  ig2_mean_out_degree <- mean(igraph::degree(ig2, mode = "out"))

  expect_lt(95, ig2_mean_out_degree)
  expect_lt(ig2_mean_out_degree, 105)


  ig3 <- sample_igraph(dfm3)
  ig3_density <- igraph::ecount(ig3) / as.numeric(n * n)

  expect_lt(0.008, ig3_density)
  expect_lt(ig3_density, 0.012)


  ### tidygraph tests ----------------------------------------------------------

  tg <- sample_tidygraph(dfm)
  tg_mean_in_degree <- mean(igraph::degree(tg, mode = "in"))

  expect_lt(9, tg_mean_in_degree)
  expect_lt(tg_mean_in_degree, 11)


  tg2 <- sample_tidygraph(dfm2)
  tg2_mean_out_degree <- mean(igraph::degree(tg2, mode = "out"))

  expect_lt(95, tg2_mean_out_degree)
  expect_lt(tg2_mean_out_degree, 105)


  tg3 <- sample_tidygraph(dfm3)
  tg3_density <- igraph::ecount(tg3) / as.numeric(n * n)

  expect_lt(0.008, tg3_density)
  expect_lt(tg3_density, 0.012)


  ### decomposition sanity check -----------------------------------------------

  expect_silent(svds(dfm))
})
