test_that("sampling from undirected factor models doesn't drop isolated nodes", {
  set.seed(27)

  latent <- sbm(
    n = 10000,
    B = diag(10),
    pi = rep(1, 10),
    expected_degree = 1,
    sort_nodes = TRUE
  )

  # when sampling edgelists, node ids for isolated nodes do not appear in the
  # edgelist. since we construct downstream network objects from the edgelists
  # this can lead to issues where sampled networks accidentally drop isolated
  # nodes

  el <- sample_edgelist(latent)

  expect_lte(
    length(unique(el$from)), # not n! and that's okay!,
    10000
  )

  set.seed(27)

  A <- sample_sparse(latent)

  expect_equal(
    dim(A),
    c(10000, 10000)
  )

  set.seed(27)

  ig <- sample_igraph(latent)

  expect_equal(
    igraph::vcount(ig),
    10000
  )

  set.seed(27)

  tbl_graph <- sample_tidygraph(latent)

  expect_equal(
    igraph::vcount(tbl_graph),
    10000
  )
})

test_that("sampling from square directed factor models doesn't drop isolated nodes", {
  set.seed(32)

  bm <- as.matrix(cbind(
    c(.3, .005, .005, .005, .005),
    c(.002, .3, .005, .005, .005),
    c(.002, .01, .3, .005, .005),
    c(.002, .01, .005, .2, .005),
    c(.002, .005, .005, .005, .2)
  ))

  pi <- c(5, 50, 20, 25, 100)

  latent <- fastRG::directed_dcsbm(
    B = bm,
    pi_in = pi,
    pi_out = pi,
    theta_in = rep(1, 200),
    theta_out = rep(1, 200),
    expected_out_degree = 3,
    allow_self_loops = FALSE,
    sort_nodes = TRUE
  )

  el <- sample_edgelist(latent)

  expect_lte(
    length(unique(c(el$from, el$to))), # not n! and that's okay!,
    200
  )

  set.seed(32)

  A <- sample_sparse(latent)

  expect_equal(
    dim(A),
    c(200, 200)
  )

  set.seed(32)

  ig <- sample_igraph(latent)

  expect_equal(
    igraph::vcount(ig),
    200
  )

  set.seed(32)

  tbl_graph <- sample_tidygraph(latent)

  expect_equal(
    igraph::vcount(tbl_graph),
    200
  )
})


test_that("sampling from rectangular directed factor models doesn't drop isolated nodes", {
  n <- 10000

  k1 <- 5
  k2 <- 3

  d <- 5000

  X <- matrix(rpois(n = n * k1, 1), nrow = n)
  S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
  Y <- matrix(rexp(n = k2 * d, 1), nrow = d)

  latent <- directed_factor_model(
    X,
    S,
    Y,
    expected_in_degree = 1
  )

  el <- sample_edgelist(latent)

  # nodes with out-degree > 0
  expect_lte(
    length(unique(el$from)), # not n! and that's okay!,
    n
  )

  # nodes with in-degree > 0
  expect_lte(
    length(unique(el$to)), # not d! and that's okay!,
    d
  )

  set.seed(32)

  A <- sample_sparse(latent)

  expect_equal(
    dim(A),
    c(n, d)
  )

  set.seed(32)

  ig <- sample_igraph(latent)

  # total node count for bipartite graph should be n + d
  expect_equal(
    igraph::vcount(ig),
    n + d
  )

  # one mode of the graph
  expect_equal(
    sum(!igraph::V(ig)$type),
    n
  )

  # the other mode of the graph
  expect_equal(
    sum(igraph::V(ig)$type),
    d
  )

  set.seed(32)

  tbl_graph <- sample_tidygraph(latent)

  # total node count for bipartite graph should be n + d
  expect_equal(
    igraph::vcount(tbl_graph),
    n + d
  )

  # one mode of the graph
  expect_equal(
    sum(!igraph::V(tbl_graph)$type),
    n
  )

  # the other mode of the graph
  expect_equal(
    sum(igraph::V(tbl_graph)$type),
    d
  )
})
