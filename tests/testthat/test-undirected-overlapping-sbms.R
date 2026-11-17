test_that("rank 1 overlapping sbms sample", {
  set.seed(27)

  n <- 10
  k <- 1
  pi <- rep(1, k) / k

  B <- matrix(0.5)

  sbm <- overlapping_sbm(n = n, pi = pi, B = B)

  expect_silent(A <- sample_sparse(sbm))
})

test_that("every node belongs to at least one block", {
  set.seed(27)

  # with these pi, about (1 - 0.1)^3 = 73% of nodes join no block at first
  osbm <- overlapping_sbm(
    n = 1000,
    B = diag(0.8, 3) + 0.1,
    pi = c(0.1, 0.1, 0.1),
    force_pure = FALSE,
    expected_density = 0.01
  )

  expect_true(all(Matrix::rowSums(osbm$Z) >= 1))
  expect_identical(osbm$X, osbm$Z)
})

test_that("nodes in no block are reassigned proportionally to pi", {
  set.seed(27)

  # block 1 is never joined directly, so every node in it was reassigned,
  # and the other two blocks are rarely joined directly
  pi <- c(0, 0.01, 0.03)

  osbm <- overlapping_sbm(
    n = 20000,
    B = diag(0.8, 3) + 0.1,
    pi = pi,
    force_pure = FALSE,
    sort_nodes = FALSE,
    expected_density = 0.01
  )

  Z <- as.matrix(osbm$Z)

  expect_true(all(rowSums(Z) >= 1))
  expect_equal(sum(Z[, 1]), 0)

  # the reassigned blocks follow pi / sum(pi) = (0, 0.25, 0.75), and
  # reassigned nodes make up the large majority of all nodes
  single <- Z[rowSums(Z) == 1, , drop = FALSE]
  props <- colMeans(single)
  expect_equal(unname(props), c(0, 0.25, 0.75), tolerance = 0.05)
})

test_that("pi must have a positive element", {
  expect_error(
    overlapping_sbm(n = 10, B = diag(2), pi = c(0, 0)),
    "must be positive"
  )
})

test_that("rank 1 overlapping sbms put every node in the block", {
  set.seed(27)

  osbm <- overlapping_sbm(n = 50, B = matrix(0.5), pi = 0.2)
  expect_true(all(osbm$Z == 1))
})

test_that("poisson_edges and allow_self_loops are passed through", {
  set.seed(27)

  B <- diag(0.8, 3) + 0.1

  osbm <- overlapping_sbm(
    n = 200,
    B = B,
    poisson_edges = FALSE,
    allow_self_loops = FALSE,
    expected_density = 0.2
  )

  expect_false(osbm$poisson_edges)
  expect_false(osbm$allow_self_loops)

  # S is back-transformed to a Poisson parameterization
  poisson_osbm <- overlapping_sbm(
    n = 200,
    B = B,
    expected_density = 0.2
  )
  expect_true(poisson_osbm$poisson_edges)
  expect_true(poisson_osbm$allow_self_loops)

  A <- sample_sparse(osbm)
  expect_equal(sum(Matrix::diag(A)), 0)
  expect_lte(max(A), 1)

  el <- sample_edgelist(osbm)
  expect_false(any(el$from == el$to))
  expect_false(anyDuplicated(el) > 0)
})

test_that("Bernoulli overlapping SBMs follow the noisy-OR edge probability", {
  set.seed(27)

  # with pi = c(1, 1) and no pure nodes, every node is in both blocks
  B <- matrix(c(0.3, 0.05, 0.05, 0.2), 2, 2)

  osbm <- overlapping_sbm(
    n = 300,
    B = B,
    pi = c(1, 1),
    poisson_edges = FALSE,
    allow_self_loops = FALSE,
    force_pure = FALSE
  )

  # every node is in both blocks, so P(edge) = 1 - prod (1 - B_ab)
  expect_true(all(osbm$Z == 1))
  p <- 1 - prod(1 - B)

  reps <- 20
  density <- mean(replicate(reps, {
    A <- sample_sparse(osbm)
    sum(A) / (300 * 299)
  }))

  expect_equal(density, p, tolerance = 0.02)
})

test_that("k and B cannot both be specified", {
  expect_error(
    overlapping_sbm(n = 10, k = 2, B = diag(0.5, 2)),
    "not both"
  )
})

test_that("print() reports whether nodes are sorted", {
  set.seed(27)

  B <- diag(0.8, 2) + 0.1

  sorted <- overlapping_sbm(n = 20, B = B, sort_nodes = TRUE)
  unsorted <- overlapping_sbm(n = 20, B = B, sort_nodes = FALSE)

  expect_output(print(sorted), "Nodes \\(n\\): 20 \\(sorted\\)")
  expect_output(print(unsorted), "Nodes \\(n\\): 20 \\(not sorted\\)")
  expect_output(print(sorted), "Undirected Overlapping Stochastic Blockmodel")
})
