library(magrittr)

test_that("undirected expected degree computed consistently", {

  # see issue 19

  set.seed(27)

  n <- 1000
  pop <- n / 2
  a <- .1
  b <- .05
  B <- matrix(c(a, b, b, a), nrow = 2)

  b_model <- sbm(
    n = n, k = 2,
    B = B,
    poisson_edges = FALSE
  )

  expect_equal(
    expected_degree(b_model), # computed
    pop * a + pop * b,       # expected "undirected edge degree",
    tolerance = 5
  )

  A <- sample_sparse(b_model)

  ### degree computation gotchas

  mean(rowSums(A))  # double counts undirected edges
  #> [1] 156.711
  mean(rowSums(triu(A)))  # right way to count undirected edges in A
  #> [1] 78.413

  expect_equal(
    mean(rowSums(triu(A))), # computed
    pop * a + pop * b,      # expected "undirected edge degree"
    tolerance = 5
  )

  model2 <- sbm(n = n, k = 2, B = B, poisson_edges = FALSE, expected_degree = 75)

  expect_equal(
    expected_degree(model2), # computed
    pop * a + pop * b,       # expected "undirected edge degree",
    tolerance = 5
  )

  A2 <- sample_sparse(model2)

  expect_equal(
    mean(rowSums(triu(A2))), # computed
    pop * a + pop * b,       # expected "undirected edge degree",
    tolerance = 5
  )

})

test_that("undirected density computed consistently", {

  # see issue 19

  set.seed(27)

  n <- 1000
  pop <- n / 2
  a <- .1
  b <- .05
  B <- matrix(c(a, b, b, a), nrow = 2)

  b_model <- sbm(
    n = n, k = 2,
    B = B,
    poisson_edges = FALSE
  )

  expect_equal(
    expected_density(b_model),              # computed
    n * (pop * a + pop * b) / choose(n, 2), # expected undirected degree density, possibly being a little sloppy about the diagonal
    tolerance = 0.05
  )

  A <- sample_sparse(b_model)

  ### density computation gotchas

  # almost correct because double counts UT and LT in num and denom,
  # but diagonal gets too much weight. slight over-estimate of density
  sum(A) / n^2

  sum(triu(A)) / choose(n, 2)  # correct density estimate

  expect_equal(
    sum(triu(A)) / choose(n, 2),            # computed
    n * (pop * a + pop * b) / choose(n, 2), # expected "undirected edge degree",
    tolerance = 0.05
  )

  model2 <- sbm(n = n, k = 2, B = B, expected_density = 0.15)

  expect_equal(
    expected_density(model2),              # computed
    n * (pop * a + pop * b) / choose(n, 2), # expected undirected degree density, possibly being a little sloppy about the diagonal
    tolerance = 0.02
  )

  A2 <- sample_sparse(model2)

  expect_equal(
    sum(triu(A2)) / choose(n, 2),            # computed
    0.15, # expected "undirected edge degree",
    tolerance = 0.05
  )

})

test_that("undirected factor model", {

  library(tidygraph)

  set.seed(7)

  n <- 1000
  k <- 5

  # don't allow self edges at all in these calculations via
  # SBM model with zero on the diagonal of B

  B <- matrix(data = 0.5, nrow = k, ncol = k)
  diag(B) <- 0

  ufm <- sbm(n = n, k = k, B = B, expected_degree = 10)
  expect_equal(expected_degree(ufm), 10)
  expect_equal(expected_density(ufm), 0.02, tolerance = 0.05) # tolerance should be relative here


  el <- sample_edgelist(ufm)
  el_mean_degree <- nrow(el) / n
  expect_lt(9, el_mean_degree)
  expect_lt(el_mean_degree, 11)

  g2 <- igraph::graph_from_data_frame(el, directed = TRUE)
  A <- igraph::as_adj(g2)

  # NOTE: see issue #19 about the following
  #
  # mean(rowSums(A))        # double counts undirected edges
  # mean(rowSums(triu(A)))  # right way to count undirected edges

  A <- sample_sparse(ufm)
  matrix_mean_degree <- mean(rowSums(triu(A)))

  expect_equal(rowSums(A), colSums(A))
  expect_lt(9, matrix_mean_degree)
  expect_lt(matrix_mean_degree, 11)


  graph <- sample_igraph(ufm)
  # igraph doubles edge counts relative to the way we want to count
  igraph_mean_degree <- mean(igraph::degree(graph)) / 2

  expect_lt(9, igraph_mean_degree)
  expect_lt(igraph_mean_degree, 11)


  tbl_graph <- sample_tidygraph(ufm)
  tbl_graph_edges <- tbl_graph %>%
    activate(edges) %>%
    as_tibble() %>%
    nrow()

  tbl_graph_mean_degree <- tbl_graph_edges / n

  expect_lt(9, tbl_graph_mean_degree)
  expect_lt(tbl_graph_mean_degree, 11)

  expect_silent(eigs_sym(ufm))
})

test_that("directed factor model", {

  set.seed(8)

  library(dplyr)

  n <- 5000
  d <- 800

  k1 <- 5
  k2 <- 3

  X <- matrix(rpois(n = n * k1, 1), nrow = n)
  Y <- matrix(rpois(n = d * k2, 1), nrow = d)
  S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1)

  dfm <- directed_factor_model(X = X, S = S, Y = Y, expected_in_degree = 10)
  expect_equal(expected_in_degree(dfm), 10)

  dfm2 <- directed_factor_model(X = X, S = S, Y = Y, expected_out_degree = 100)
  expect_equal(expected_out_degree(dfm2), 100)

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

  el3_density <- nrow(el3) / as.numeric(n * d)

  expect_lt(0.08, el3_density)
  expect_lt(el3_density, 0.12)


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

  expect_lt(0.08, A3_density)
  expect_lt(A3_density, 0.12)

  ### igraph tests --------------------------------------------------------------

  ig <- sample_igraph(dfm)
  A_ig <- igraph::as_incidence_matrix(ig, sparse = TRUE, names = FALSE)
  ig_mean_in_degree <- mean(colSums(A_ig))

  expect_lt(9, ig_mean_in_degree)
  expect_lt(ig_mean_in_degree, 11)


  ig2 <- sample_igraph(dfm2)
  A2_ig <- igraph::as_incidence_matrix(ig2, sparse = TRUE, names = FALSE)
  ig2_mean_out_degree <- mean(rowSums(A2_ig))


  expect_lt(95, ig2_mean_out_degree)
  expect_lt(ig2_mean_out_degree, 105)


  ig3 <- sample_igraph(dfm3)
  ig3_density <- igraph::ecount(ig3) / as.numeric(n * d)

  expect_lt(0.08, ig3_density)
  expect_lt(ig3_density, 0.12)


  ### tidygraph tests ----------------------------------------------------------

  tg <- sample_tidygraph(dfm)
  A_tg <- igraph::as_incidence_matrix(tg, sparse = TRUE, names = FALSE)
  tg_mean_in_degree <- mean(colSums(A_tg))

  expect_lt(9, tg_mean_in_degree)
  expect_lt(tg_mean_in_degree, 11)


  tg2 <- sample_tidygraph(dfm2)
  A2_tg <- igraph::as_incidence_matrix(tg2, sparse = TRUE, names = FALSE)
  tg2_mean_out_degree <- mean(rowSums(A2_tg))


  expect_lt(95, tg2_mean_out_degree)
  expect_lt(tg2_mean_out_degree, 105)


  tg3 <- sample_tidygraph(dfm3)
  tg3_density <- igraph::ecount(tg3) / as.numeric(n * d)

  expect_lt(0.08, tg3_density)
  expect_lt(tg3_density, 0.12)


  ### decomposition sanity check -----------------------------------------------

  expect_silent(svds(dfm))
})
