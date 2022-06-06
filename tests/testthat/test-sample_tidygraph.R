test_that("doesn't drop unconnected nodes", {

  set.seed(27)

  latent <- sbm(
    n = 10000,
    B = diag(10),
    pi = rep(1, 10),
    expected_degree = 1,
    sort_nodes = TRUE
  )

  graph <- sample_tidygraph(latent)

  expect_equal(
    igraph::vcount(graph),
    10000
  )

})
