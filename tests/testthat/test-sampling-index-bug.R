test_that("degenerate test case 1 succeeds", {
  set.seed(7)

  n <- 10
  k <- 10
  pi <- rep(1, k) / k

  B <- diag(rep(0.5, k))

  sbm <- sbm(n = n, pi = pi, B = B)

  expect_silent(sample_sparse(sbm))
})

test_that("degenerate test case 2 succeeds", {
  set.seed(7)

  n <- 10
  k <- 10
  pi <- rep(1, k) / k

  B <- diag(rep(0.5, k))

  dsbm <- directed_dcsbm(pi_in = pi, pi_out = pi, B = B, theta_in = rep(1, n), theta_out = rep(1, n), sort_nodes = TRUE)

  expect_silent(sample_sparse(dsbm))
})

test_that("degenerate test case 3 succeeds", {
  set.seed(9)

  n <- 10
  d <- 5
  k_in <- 10
  k_out <- 5

  pi_in <- rep(1, k_in)
  pi_out <- rep(1, k_out)

  B <- matrix(0, nrow = k_in, ncol = k_out)
  diag(B) <- 0.5

  dsbm <- directed_dcsbm(pi_in = pi_in, pi_out = pi_out, B = B, theta_in = rep(1, n), theta_out = rep(1, n), sort_nodes = TRUE)

  expect_silent(sample_sparse(dsbm))
})
