library(fastRG)

set.seed(27)

test_that("doesn't fail // sufficient dummy columns in SBMs", {

  n <- 10
  k <- 10
  pi <- rep(1, k) / k

  B <- diag(rep(0.5, k))

  sbm <- sbm(n = n, pi = pi, B = B)

  expect_silent(sample_sparse(sbm))
})
