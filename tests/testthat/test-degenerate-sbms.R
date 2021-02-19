test_that("degenerate SBMs", {

  set.seed(5)

  n <- 10
  k <- 10
  pi <- rep(1, k) / k

  B <- diag(rep(0.5, k))

  sbm <- sbm(n = n, pi = pi, B = B)

  expect_silent(sample_sparse(sbm))
})
