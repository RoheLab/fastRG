# shared helpers for property-based tests. these check invariants that
# must hold for any model, rather than pinning particular values

# column of the (single) nonzero entry in each row of a block membership
# matrix, with 0 for rows that are entirely zero
membership_column <- function(X) {
  X <- as.matrix(X)
  out <- max.col(X != 0, ties.method = "first")
  out[rowSums(X != 0) == 0] <- 0L
  out
}

# dense P = E[A | X, S] for any factor model
dense_expectation <- function(model) {
  as.matrix(expectation(model))
}

# z-score of an observed total of `reps` independent Poisson(mu) draws
poisson_total_z <- function(total, mu, reps = 1) {
  (total - reps * mu) / sqrt(reps * mu)
}

# Bonferroni-corrected two-sided normal critical value for `m` simultaneous
# z-tests at overall level `alpha`
bonferroni_z <- function(m, alpha = 1e-3) {
  stats::qnorm(1 - alpha / (2 * m))
}

# check the internal consistency of a (undirected) dcsbm-like model: z, X,
# theta and pi all describe the same nodes, in the same order
expect_dcsbm_consistent <- function(m) {
  X <- as.matrix(m$X)
  n <- m$n
  k <- m$k

  expect_equal(dim(X), c(n, k))
  expect_s3_class(m$z, "factor")
  expect_length(m$z, n)
  expect_equal(nlevels(m$z), k)
  expect_length(m$theta, n)
  expect_length(m$pi, k)
  expect_equal(sum(m$pi), 1)
  expect_true(all(m$pi >= 0))

  # each row of X has exactly one nonzero entry, in column z[i], and it is
  # theta[i]
  expect_true(all(rowSums(X != 0) <= 1))
  nonzero <- m$theta > 0
  expect_equal(membership_column(X)[nonzero], as.integer(m$z)[nonzero])
  expect_equal(unname(rowSums(X)), as.numeric(m$theta))

  # S is symmetric and the right size
  expect_equal(dim(m$S), c(k, k))
  expect_true(isSymmetric(as.matrix(m$S)))

  invisible(m)
}

# check the internal consistency of a directed_dcsbm
expect_directed_dcsbm_consistent <- function(m) {
  X <- as.matrix(m$X)
  Y <- as.matrix(m$Y)

  expect_equal(dim(X), c(m$n, m$k1))
  expect_equal(dim(Y), c(m$d, m$k2))
  expect_equal(dim(m$S), c(m$k1, m$k2))

  expect_length(m$z_out, m$n)
  expect_length(m$z_in, m$d)
  expect_equal(nlevels(m$z_out), m$k1)
  expect_equal(nlevels(m$z_in), m$k2)
  expect_length(m$theta_out, m$n)
  expect_length(m$theta_in, m$d)
  expect_length(m$pi_out, m$k1)
  expect_length(m$pi_in, m$k2)
  expect_equal(sum(m$pi_out), 1)
  expect_equal(sum(m$pi_in), 1)

  expect_equal(membership_column(X), as.integer(m$z_out))
  expect_equal(membership_column(Y), as.integer(m$z_in))
  expect_equal(unname(rowSums(X)), as.numeric(m$theta_out))
  expect_equal(unname(rowSums(Y)), as.numeric(m$theta_in))

  invisible(m)
}

# all four samplers, run from the same seed. they all consume the RNG the
# same way, so they describe the same graph
sample_all <- function(model, seed, multigraph = TRUE) {
  set.seed(seed)
  el <- sample_edgelist(model, multigraph = multigraph)
  set.seed(seed)
  A <- sample_sparse(model)
  set.seed(seed)
  ig <- sample_igraph(model, multigraph = multigraph)
  set.seed(seed)
  tg <- sample_tidygraph(model, multigraph = multigraph)
  list(el = el, A = A, ig = ig, tg = tg)
}
