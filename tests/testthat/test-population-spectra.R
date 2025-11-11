# sin Theta dist in Frob norm
sin_theta_distance <- function(u, v) {
  s <- svd(crossprod(u, v))
  ncol(u) - sum(s$d^2)
}

test_that("tested eigs_sym and svds for population spectra on sbm with repeated eigenvalues", {
  set.seed(27)

  n <- 100
  k <- 5
  B <- diag(rep(0.5, k))

  sbm <- sbm(n = n, B = B, allow_self_loops = FALSE)

  EA <- expectation(sbm)
  EA_manual <- sbm$X %*% tcrossprod(sbm$S, sbm$X)

  expect_equal(
    expected_degree(sbm),
    10
  )

  expect_equal(EA_manual, EA)

  eig_manual <- unclass(eigen(EA))
  eig_manual$values <- eig_manual$values[1:k]
  eig_manual$vectors <- eig_manual$vectors[, 1:k]
  eig <- eigs_sym(sbm)
  s <- svds(sbm)

  expect_equal(eig$values, eig_manual$values)
  expect_equal(s$d, eig_manual$values)

  expect_equal(0, sin_theta_distance(eig$vectors, eig_manual$vectors))
  expect_equal(0, sin_theta_distance(s$u, eig_manual$vectors))
  expect_equal(0, sin_theta_distance(s$v, eig_manual$vectors))

  el <- sample_edgelist(sbm)
  A <- sample_sparse(sbm)
  s_obs <- svds(A, k)

  expect_true(
    all(s_obs$d >= s$d - 2 * log(n) & s_obs$d <= s$d + 2 * log(n))
  )
})

test_that("tested eigs_sym and svds for population spectra on sbm with distinct eigenvalues", {
  set.seed(27)

  n <- 100
  k <- 10
  B <- matrix(0.03, nrow = k, ncol = k)
  diag(B) <- 0.5

  sbm <- dcsbm(B = B, theta = rexp(n) + 1)

  EA <- expectation(sbm)
  EA_manual <- sbm$X %*% tcrossprod(sbm$S, sbm$X)

  expect_equal(EA_manual, EA)

  eig_manual <- unclass(eigen(EA))
  eig_manual$values <- eig_manual$values[1:k]
  eig_manual$vectors <- eig_manual$vectors[, 1:k]
  eig <- eigs_sym(sbm)
  s <- svds(sbm)

  expect_equal(eig$values, eig_manual$values)
  expect_equal(s$d, eig_manual$values)

  expect_equal(0, sin_theta_distance(eig$vectors, eig_manual$vectors))
  expect_equal(0, sin_theta_distance(s$u, eig_manual$vectors))
  expect_equal(0, sin_theta_distance(s$v, eig_manual$vectors))

  A <- sample_sparse(sbm)
  s_obs <- irlba::irlba(A, k)

  expect_true(
    all(s_obs$d >= s$d - 2 * log(n) & s_obs$d <= s$d + 2 * log(n))
  )
})
