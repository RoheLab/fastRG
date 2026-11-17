# see issue #46: undirected samplers used to draw self-loops at half the
# intended rate, a side effect of passing S / 2 to the sampler to fix #19.
# these tests check that E[A | X, S] = X S X' holds on *and* off the diagonal,
# so that neither the diagonal undersampling (#46) nor the off-diagonal
# oversampling (#19, #31) creeps back in

# entrywise moments of `reps` sampled adjacency matrices, as dense matrices
adjacency_moments <- function(model, reps) {
  A_sum <- matrix(0, nrow = model$n, ncol = model$n)
  A_sq_sum <- A_sum
  A_zero <- A_sum
  A_max <- 0

  for (r in seq_len(reps)) {
    A <- as.matrix(sample_sparse(model))
    A_sum <- A_sum + A
    A_sq_sum <- A_sq_sum + A^2
    A_zero <- A_zero + (A == 0)
    A_max <- max(A_max, max(A))
  }

  A_mean <- A_sum / reps

  list(
    mean = A_mean,
    var = A_sq_sum / reps - A_mean^2,
    zero = A_zero / reps,
    max = A_max
  )
}

# bernoulli sample_sparse() output is a pattern matrix, which silently
# collapses duplicate edges, so also look at the edgelists themselves. returns
# the mean number of self-loop rows per node and whether any edgelist had
# duplicate rows
edgelist_self_loops <- function(model, reps) {
  loop_counts <- numeric(model$n)
  any_duplicates <- FALSE

  for (r in seq_len(reps)) {
    el <- sample_edgelist(model)
    loops <- el$from[el$from == el$to]
    loop_counts <- loop_counts + tabulate(loops, nbins = model$n)
    any_duplicates <- any_duplicates || anyDuplicated(el) > 0
  }

  list(mean = loop_counts / reps, any_duplicates = any_duplicates)
}

off_diagonal_mean <- function(A) {
  (sum(A) - sum(diag(A))) / (nrow(A) * (nrow(A) - 1))
}

test_that("undirected poisson graphs sample A_ii ~ Poisson((X S X')_ii)", {
  set.seed(46)

  n <- 100
  lambda <- 0.4

  # E[A] = lambda in every entry, including the diagonal
  ufm <- undirected_factor_model(X = matrix(1, n, 1), S = matrix(lambda, 1, 1))

  moments <- adjacency_moments(ufm, reps = 500)

  # 50,000 Poisson(0.4) draws on the diagonal, so sd of the mean is ~0.003.
  # undersampling gives 0.2 and oversampling gives 0.8
  expect_equal(mean(diag(moments$mean)), lambda, tolerance = 0.025)
  expect_equal(off_diagonal_mean(moments$mean), lambda, tolerance = 0.01)

  # the mean alone can't tell Poisson(lambda) apart from, e.g., doubling a
  # Poisson(lambda / 2), which has variance 2 * lambda and
  # P(0) = exp(-lambda / 2), so also check the variance and P(A_ii = 0)
  expect_equal(mean(diag(moments$var)), lambda, tolerance = 0.05)
  expect_equal(mean(diag(moments$zero)), exp(-lambda), tolerance = 0.01)
})

test_that("undirected bernoulli graphs sample A_ii ~ Bernoulli(p_ii)", {
  set.seed(47)

  n <- 100
  p <- 0.4

  ufm <- undirected_factor_model(
    X = matrix(1, n, 1),
    S = matrix(p, 1, 1),
    poisson_edges = FALSE
  )

  moments <- adjacency_moments(ufm, reps = 500)

  # undersampling gives 1 - sqrt(1 - p) ~= 0.225 here
  expect_equal(mean(diag(moments$mean)), p, tolerance = 0.025)
  expect_equal(off_diagonal_mean(moments$mean), p, tolerance = 0.01)

  # no duplicate edges, self-loops included, in any edgelist or igraph. this
  # fails if duplicates are removed before the extra self-loops are added
  loops <- edgelist_self_loops(ufm, reps = 500)
  expect_false(loops$any_duplicates)
  expect_equal(mean(loops$mean), p, tolerance = 0.025)
  expect_false(igraph::any_multiple(sample_igraph(ufm)))
})

for (poisson_edges in c(TRUE, FALSE)) {
  edge_type <- if (poisson_edges) "poisson" else "bernoulli"

  test_that(
    paste("sbm diagonal matches B within each block for", edge_type, "edges"),
    {
      set.seed(53)

      # sbm() uses indicator X with theta = 1, so in both edge models
      # E[A_ij] = B[z_i, z_j] exactly, including for i = j
      B <- matrix(
        c(
          0.6, 0.2, 0.05,
          0.2, 0.3, 0.1,
          0.05, 0.1, 0.1
        ),
        nrow = 3
      )

      model <- sbm(
        n = 90,
        B = B,
        block_sizes = c(30, 30, 30),
        poisson_edges = poisson_edges
      )

      z <- max.col(as.matrix(model$X))
      moments <- adjacency_moments(model, reps = 1000)

      # 30,000 diagonal draws per block, so sd of each block mean is < 0.005.
      # undersampling is off by at least 0.049 (block 3, bernoulli)
      diag_by_block <- tapply(diag(moments$mean), z, mean)
      expect_lt(max(abs(diag_by_block - diag(B))), 0.015)

      off_diag <- row(moments$mean) != col(moments$mean)
      for (a in 1:3) {
        for (b in 1:3) {
          in_blocks <- off_diag & outer(z == a, z == b)
          expect_lt(abs(mean(moments$mean[in_blocks]) - B[a, b]), 0.01)
        }
      }

      # self-loop rates in the edgelists match too, not just in sample_sparse()
      loops <- edgelist_self_loops(model, reps = 1000)
      loops_by_block <- tapply(loops$mean, z, mean)
      expect_lt(max(abs(loops_by_block - diag(B))), 0.015)

      if (poisson_edges) {
        var_by_block <- tapply(diag(moments$var), z, mean)
        expect_lt(max(abs(var_by_block - diag(B))), 0.03)
      } else {
        expect_false(loops$any_duplicates)
      }
    }
  )
}

test_that("sampled diagonal matches diag(X S X') for heterogeneous factors", {
  set.seed(48)

  n <- 60
  k <- 3

  X <- matrix(rexp(n * k), nrow = n)
  S <- matrix(c(0.6, 0.1, 0.2, 0.1, 0.5, 0.1, 0.2, 0.1, 0.4), nrow = k)

  ufm <- undirected_factor_model(X, S)
  EA <- as.matrix(expectation(ufm))

  A_bar <- adjacency_moments(ufm, reps = 1000)$mean

  expect_equal(sum(diag(A_bar)) / sum(diag(EA)), 1, tolerance = 0.03)
  expect_equal(diag(A_bar), diag(EA), tolerance = 0.2)

  off_diag <- row(EA) != col(EA)
  expect_equal(sum(A_bar[off_diag]) / sum(EA[off_diag]), 1, tolerance = 0.02)

  # with the diagonal fixed, sampled totals match expected_edges() and
  # expected_degree() exactly in expectation (#19)
  expect_equal(sum(A_bar) / expected_edges(ufm), 1, tolerance = 0.02)
  expect_equal(mean(rowSums(A_bar)), expected_degree(ufm), tolerance = 0.02)

  # average of sampled graphs has the population spectrum, rather than
  # being off by a factor of two (#31)
  pop_values <- eigs_sym(ufm)$values
  sample_values <- eigen(A_bar, symmetric = TRUE)$values[1:k]
  expect_equal(sample_values, pop_values, tolerance = 0.03)
})

test_that("sampling does not modify the model", {
  set.seed(49)

  B <- matrix(0.1, nrow = 3, ncol = 3)
  diag(B) <- 0.5

  model <- sbm(n = 200, B = B, expected_degree = 10)
  S_before <- model$S

  sample_sparse(model)
  sample_edgelist(model)

  expect_equal(model$S, S_before)
  expect_equal(expected_degree(model), 10)
})

for (poisson_edges in c(TRUE, FALSE)) {
  edge_type <- if (poisson_edges) "poisson" else "bernoulli"

  test_that(
    paste("allow_self_loops = FALSE removes all self-loops for", edge_type, "edges"),
    {
      set.seed(50)

      n <- 100
      lambda <- 0.4

      ufm <- undirected_factor_model(
        X = matrix(1, n, 1),
        S = matrix(lambda, 1, 1),
        poisson_edges = poisson_edges,
        allow_self_loops = FALSE
      )

      el <- sample_edgelist(ufm)
      expect_false(any(el$from == el$to))

      A_bar <- adjacency_moments(ufm, reps = 200)$mean
      expect_equal(sum(diag(A_bar)), 0)
      expect_equal(off_diagonal_mean(A_bar), lambda, tolerance = 0.02)
    }
  )
}

test_that("undirected edgelists keep from <= to with self-loops added", {
  set.seed(51)

  ufm <- undirected_factor_model(X = matrix(1, 50, 1), S = matrix(0.5, 1, 1))
  el <- sample_edgelist(ufm)

  expect_true(any(el$from == el$to))
  expect_true(all(el$from <= el$to))
})

test_that("undirected models with no edges still sample", {
  ufm <- undirected_factor_model(X = matrix(1, 10, 1), S = matrix(0, 1, 1))

  expect_equal(nrow(sample_edgelist(ufm)), 0)
  expect_equal(sum(sample_sparse(ufm)), 0)
  expect_silent(sample_igraph(ufm))

  bern <- undirected_factor_model(
    X = matrix(1, 10, 1),
    S = matrix(0, 1, 1),
    poisson_edges = FALSE,
    allow_self_loops = FALSE
  )

  expect_equal(nrow(sample_edgelist(bern)), 0)
  expect_equal(sum(sample_sparse(bern)), 0)
})

test_that("directed graphs with square expectation sample the diagonal correctly", {
  set.seed(52)

  n <- 100
  lambda <- 0.4

  dfm <- directed_factor_model(
    X = matrix(1, n, 1),
    S = matrix(lambda, 1, 1),
    Y = matrix(1, n, 1)
  )

  A_bar <- adjacency_moments(dfm, reps = 500)$mean

  expect_equal(mean(diag(A_bar)), lambda, tolerance = 0.025)
  expect_equal(off_diagonal_mean(A_bar), lambda, tolerance = 0.01)
})
