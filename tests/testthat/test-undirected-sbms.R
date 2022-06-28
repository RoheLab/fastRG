test_that("SBMs don't drop isolated nodes", {

  set.seed(27)

  n <- 10
  k <- 10
  pi <- rep(1, k) / k

  B <- diag(rep(0.5, k))

  sbm <- sbm(n = n, pi = pi, B = B)

  A <- sample_sparse(sbm)

  expect_equal(ncol(A), 10)
  expect_equal(nrow(A), 10)
})

test_that("rank 1 sbms sample", {

  set.seed(27)

  n <- 10
  k <- 1
  pi <- rep(1, k) / k

  B <- matrix(0.5)

  sbm <- sbm(n = n, pi = pi, B = B)

  expect_silent(A <- sample_sparse(sbm))
})


