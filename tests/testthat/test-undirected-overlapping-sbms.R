test_that("rank 1 overlapping sbms sample", {

  set.seed(27)

  n <- 10
  k <- 1
  pi <- rep(1, k) / k
  theta <- rep(1, n)

  B <- matrix(0.5)

  sbm <- overlapping_sbm(n = n, pi = pi, B = B, theta = theta)

  expect_silent(A <- sample_sparse(sbm))
})


