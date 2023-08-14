test_that("doesn't drop isolated nodes (#35)", {

  set.seed(32)

  bm <- as.matrix(cbind(
    c(.3, .005, .005, .005, .005),
    c(.002, .3, .005, .005, .005),
    c(.002, .01, .3, .005, .005),
    c(.002, .01, .005, .2, .005),
    c(.002, .005, .005, .005, .2)
  ))

  pi <- c(5, 50, 20, 25, 100)

  sbm <- fastRG::directed_dcsbm(
    B = bm,
    pi_in = pi,
    pi_out = pi,
    theta_in = rep(1, 200),
    theta_out = rep(1, 200),
    expected_out_degree = 3,
    allow_self_loops = FALSE,
    sort_nodes = TRUE
  )

  net <- fastRG::sample_igraph(sbm)

  expect_equal(
    igraph::vcount(net),
    200
  )

})
