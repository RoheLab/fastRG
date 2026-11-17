# Monte Carlo checks that sampled graphs have the expectation the model
# claims, entry by entry.
#
# Write P = expectation(model). Conditional on the latent factors, the
# entries of A (upper triangle and diagonal for undirected graphs) are
# independent, with
#
#   Poisson graphs:    A[i, j] ~ Poisson(P[i, j])
#   Bernoulli graphs:  A[i, j] ~ Bernoulli(1 - exp(-P[i, j]))
#
# because Bernoulli graphs are sampled as Poisson graphs with rate P and
# then thresholded. So over `reps` independent replicates the entrywise
# totals T[i, j] are exactly Poisson(reps P[i, j]) or
# Binomial(reps, 1 - exp(-P[i, j])), and the Pearson statistic
#
#   sum over cells (T - E[T])^2 / Var[T]
#
# is approximately chi-squared with one degree of freedom per cell, as long
# as E[T] is not small. The models below keep every reps * P[i, j] above
# ~50, so the chi-squared approximation is accurate, and we test at level
# 1e-4. Cells where the model says A[i, j] = 0 (no self-loops) must be
# exactly zero in every replicate.

# sum of `reps` sampled adjacency matrices, as a dense matrix
sum_of_samples <- function(model, reps) {
  total <- matrix(0, model$n, if (is.null(model$d)) model$n else model$d)
  for (r in seq_len(reps)) {
    total <- total + as.matrix(sample_sparse(model))
  }
  total
}

# the cells of A that are sampled independently
independent_cells <- function(model) {
  P <- dense_expectation(model)
  if (inherits(model, "undirected_factor_model")) {
    upper.tri(P, diag = TRUE)
  } else {
    matrix(TRUE, nrow(P), ncol(P))
  }
}

expect_entrywise_expectation <- function(model, reps, label, level = 1e-4) {
  P <- dense_expectation(model)
  total <- sum_of_samples(model, reps)
  cells <- independent_cells(model)

  square <- nrow(P) == ncol(P)
  zero <- matrix(FALSE, nrow(P), ncol(P))
  if (!model$allow_self_loops && square) {
    diag(zero) <- TRUE
  }

  # structural zeros are exactly zero
  expect_true(all(total[zero] == 0), label = label)

  test <- cells & !zero

  if (model$poisson_edges) {
    mean <- reps * P[test]
    var <- mean
  } else {
    p <- 1 - exp(-P[test])
    mean <- reps * p
    var <- reps * p * (1 - p)
  }

  stat <- sum((total[test] - mean)^2 / var)
  df <- sum(test)

  expect_lt(stat, stats::qchisq(1 - level, df), label = paste(label, "chi-squared"))

  # and the statistic is not suspiciously small either, which would mean
  # that the variance is wrong
  expect_gt(stat, stats::qchisq(level, df), label = paste(label, "chi-squared"))
}

test_that("undirected Poisson graphs have E[A] = expectation(model)", {
  skip_on_cran()

  set.seed(50)
  n <- 10
  X <- matrix(runif(n * 2, 0.5, 1.5), n)
  S <- matrix(c(0.8, 0.3, 0.1, 0.6), 2) # asymmetric, gets symmetrized

  for (loops in c(TRUE, FALSE)) {
    model <- undirected_factor_model(X, S, allow_self_loops = loops)
    expect_entrywise_expectation(model, reps = 150, label = paste("loops", loops))
  }
})

test_that("undirected Bernoulli graphs have E[A] = 1 - exp(-expectation(model))", {
  skip_on_cran()

  set.seed(51)
  n <- 10

  for (loops in c(TRUE, FALSE)) {
    # a general factor model: X has unequal weights, so the edge
    # probabilities are not simply X S X'
    model <- undirected_factor_model(
      matrix(runif(n * 2, 0.5, 1.5), n),
      matrix(c(0.5, 0.2, 0.2, 0.4), 2),
      poisson_edges = FALSE,
      allow_self_loops = loops
    )
    expect_entrywise_expectation(model, reps = 400, label = paste("ufm loops", loops))

    # for blockmodels, the edge probabilities are exactly B
    B <- matrix(c(0.7, 0.2, 0.2, 0.4), 2)
    m <- sbm(n = n, B = B, block_sizes = c(4, 6), poisson_edges = FALSE, allow_self_loops = loops)
    z <- as.integer(m$z)
    expect_equal(1 - exp(-dense_expectation(m)), B[z, z], ignore_attr = TRUE)
    expect_entrywise_expectation(m, reps = 400, label = paste("sbm loops", loops))
  }
})

test_that("dense Bernoulli graphs near density one have the right density", {
  # p close to one stresses the -log(1 - p) transformation, and requires
  # many repeated Poisson edges to collapse into one
  set.seed(52)
  n <- 40
  p <- 0.97
  reps <- 20

  for (loops in c(TRUE, FALSE)) {
    model <- erdos_renyi(n = n, p = p, poisson_edges = FALSE, allow_self_loops = loops)
    pairs <- choose(n, 2) + if (loops) n else 0

    edges <- vapply(seq_len(reps), function(r) nrow(sample_edgelist(model)), numeric(1))

    expect_true(all(edges <= pairs))
    # the total is exactly Binomial(reps * pairs, p)
    z <- (sum(edges) - reps * pairs * p) / sqrt(reps * pairs * p * (1 - p))
    expect_lt(abs(z), 4.5)
  }
})

test_that("directed Poisson and Bernoulli graphs have the right expectation", {
  skip_on_cran()

  set.seed(53)
  n <- 8

  for (d in c(8, 11)) {
    for (poisson in c(TRUE, FALSE)) {
      # self-loops are only meaningful when rows and columns index the
      # same nodes
      for (loops in if (d == n) c(TRUE, FALSE) else TRUE) {
        model <- directed_factor_model(
          X = matrix(runif(n * 2, 0.5, 1.5), n),
          S = matrix(c(0.6, 0.1, 0.3, 0.2, 0.4, 0.5), 2),
          Y = matrix(runif(d * 3, 0.5, 1.5), d),
          poisson_edges = poisson,
          allow_self_loops = loops
        )
        label <- paste("d", d, "poisson", poisson, "loops", loops)
        expect_entrywise_expectation(model, reps = 200, label = label)
      }
    }
  }
})

test_that("blockmodel constructors have E[A] = expectation(model)", {
  skip_on_cran()

  set.seed(54)
  n <- 9

  models <- list(
    dcsbm_unsorted = dcsbm(
      theta = runif(n, 1, 2),
      B = matrix(c(0.5, 0.2, 0.2, 0.3), 2),
      pi = c(0.3, 0.7),
      sort_nodes = FALSE
    ),
    mmsbm = mmsbm(theta = runif(n, 1, 2), B = matrix(c(0.6, 0.2, 0.2, 0.5), 2), alpha = c(1, 2)),
    overlapping = overlapping_sbm(n = n, B = matrix(c(0.5, 0.2, 0.2, 0.4), 2), poisson_edges = FALSE),
    directed_dcsbm = directed_dcsbm(
      theta_out = runif(n, 1, 2),
      theta_in = runif(n, 1, 2),
      B = matrix(c(0.5, 0.2, 0.1, 0.4), 2),
      force_identifiability = FALSE,
      allow_self_loops = FALSE
    ),
    directed_pp = directed_planted_partition(
      n = n,
      k = 2,
      within_block = 0.6,
      between_block = 0.2,
      poisson_edges = FALSE
    )
  )

  for (name in names(models)) {
    expect_entrywise_expectation(models[[name]], reps = 250, label = name)
  }
})

test_that("expected_degrees() matches the mean sampled degree of every node", {
  skip_on_cran()

  # for Poisson graphs with self-loops, rowSums(A)[i] is a sum of
  # independent Poisson entries, so its total over reps is exactly
  # Poisson(reps * expected_degrees(model)[i])
  set.seed(55)
  n <- 80
  reps <- 100

  models <- list(
    dcsbm = dcsbm(theta = rexp(n) + 0.2, B = diag(3) * 0.5 + 0.1, expected_degree = 8),
    chung_lu = chung_lu(theta = rexp(n) + 0.2, expected_degree = 8),
    mmsbm = mmsbm(theta = rexp(n) + 0.2, B = diag(3) * 0.5 + 0.1, alpha = rep(1, 3), expected_degree = 8)
  )

  for (name in names(models)) {
    model <- models[[name]]
    expected <- expected_degrees(model)
    expect_equal(mean(expected), expected_degree(model))

    total <- numeric(n)
    for (r in seq_len(reps)) {
      total <- total + rowSums(sample_sparse(model))
    }

    z <- poisson_total_z(total, expected, reps)
    expect_lt(max(abs(z)), bonferroni_z(n), label = name)
  }
})

test_that("directed in- and out-degrees match rowSums and colSums of the expectation", {
  skip_on_cran()

  set.seed(56)
  n <- 40
  d <- 25
  reps <- 100

  model <- directed_factor_model(
    X = matrix(rexp(n * 2), n),
    S = matrix(runif(4), 2),
    Y = matrix(rexp(d * 2), d),
    expected_out_degree = 6
  )
  P <- dense_expectation(model)

  out_total <- numeric(n)
  in_total <- numeric(d)
  for (r in seq_len(reps)) {
    A <- sample_sparse(model)
    out_total <- out_total + rowSums(A)
    in_total <- in_total + colSums(A)
  }

  expect_lt(
    max(abs(poisson_total_z(out_total, rowSums(P), reps))),
    bonferroni_z(n)
  )
  expect_lt(
    max(abs(poisson_total_z(in_total, colSums(P), reps))),
    bonferroni_z(d)
  )
})
