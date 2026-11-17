test_that("directed_planted_partition() builds a constant-diagonal B", {
  set.seed(27)

  pp <- directed_planted_partition(
    n = 50,
    k = 3,
    within_block = 0.4,
    between_block = 0.1
  )

  expect_s3_class(
    pp,
    c("directed_planted_partition", "directed_dcsbm", "directed_factor_model"),
    exact = TRUE
  )

  expected_B <- matrix(0.1, 3, 3)
  diag(expected_B) <- 0.4

  expect_equal(unname(as.matrix(pp$S)), expected_B)
  expect_equal(pp$within_block, 0.4)
  expect_equal(pp$between_block, 0.1)

  expect_equal(pp$n, 50)
  expect_equal(pp$k1, 3)
  expect_equal(pp$k2, 3)

  # no degree heterogeneity
  expect_equal(pp$theta_out, rep(1, 50))
  expect_equal(pp$theta_in, rep(1, 50))
  expect_equal(pp$pi_out, rep(1 / 3, 3))
  expect_equal(pp$pi_in, rep(1 / 3, 3))

  # expected adjacency matrix is exactly B[z_out[i], z_in[j]]
  EA <- as.matrix(pp$X %*% pp$S %*% t(pp$Y))
  expect_equal(
    unname(EA),
    expected_B[as.integer(pp$z_out), as.integer(pp$z_in)]
  )
})

test_that("a and b parameterization matches within_block / between_block", {
  set.seed(27)
  from_ab <- directed_planted_partition(n = 100, k = 4, a = 20, b = 5)

  set.seed(27)
  from_probs <- directed_planted_partition(
    n = 100,
    k = 4,
    within_block = 0.2,
    between_block = 0.05
  )

  expect_equal(from_ab, from_probs)
  expect_equal(from_ab$within_block, 0.2)
  expect_equal(from_ab$between_block, 0.05)
})

test_that("in and out communities are independent by default", {
  set.seed(27)

  pp <- directed_planted_partition(
    n = 200,
    k = 3,
    within_block = 0.3,
    between_block = 0.05
  )

  expect_false(pp$shared_blocks)
  expect_false(identical(as.integer(pp$z_out), as.integer(pp$z_in)))
  expect_false(is.unsorted(pp$z_out))
  expect_false(is.unsorted(pp$z_in))
})

test_that("shared_blocks gives each node a single community", {
  for (sort_nodes in c(TRUE, FALSE)) {
    set.seed(27)

    pp <- directed_planted_partition(
      n = 60,
      k = 3,
      within_block = 0.3,
      between_block = 0.05,
      shared_blocks = TRUE,
      sort_nodes = sort_nodes
    )

    expect_true(pp$shared_blocks)
    expect_identical(pp$z_out, pp$z_in)

    # with shared blocks, E[A] is block constant with B on the diagonal
    EA <- as.matrix(pp$X %*% pp$S %*% t(pp$Y))
    same_block <- outer(pp$z_out, pp$z_in, "==")
    expect_true(all(EA[same_block] == 0.3))
    expect_true(all(EA[!same_block] == 0.05))
  }
})

test_that("block_sizes and pi apply to both outgoing and incoming blocks", {
  set.seed(27)

  sized <- directed_planted_partition(
    n = 90,
    k = 3,
    within_block = 0.3,
    between_block = 0.05,
    block_sizes = c(40, 20, 30)
  )

  expect_equal(as.vector(table(sized$z_out)), c(20, 30, 40))
  expect_equal(as.vector(table(sized$z_in)), c(20, 30, 40))

  set.seed(27)

  shared_sized <- directed_planted_partition(
    n = 90,
    k = 3,
    within_block = 0.3,
    between_block = 0.05,
    block_sizes = c(40, 20, 30),
    shared_blocks = TRUE,
    sort_nodes = FALSE
  )

  expect_identical(shared_sized$z_out, shared_sized$z_in)
  expect_equal(as.vector(table(shared_sized$z_out)), c(40, 20, 30))

  set.seed(27)

  weighted <- directed_planted_partition(
    n = 90,
    k = 3,
    within_block = 0.3,
    between_block = 0.05,
    pi = c(3, 1, 2)
  )

  expect_equal(weighted$pi_out, c(1, 2, 3) / 6)
  expect_equal(weighted$pi_in, c(1, 2, 3) / 6)
})

test_that("expected degree scaling and edge distributions pass through", {
  set.seed(27)

  pp <- directed_planted_partition(
    n = 200,
    k = 4,
    within_block = 0.5,
    between_block = 0.1,
    expected_out_degree = 5
  )

  expect_equal(expected_out_degree(pp), 5)

  # rescaling keeps B a planted partition
  S <- as.matrix(pp$S)
  expect_equal(length(unique(diag(S))), 1)
  expect_equal(length(unique(S[upper.tri(S) | lower.tri(S)])), 1)
  expect_equal(S[1, 1] / S[1, 2], 5)

  bernoulli <- directed_planted_partition(
    n = 100,
    k = 2,
    a = 30,
    b = 5,
    poisson_edges = FALSE,
    allow_self_loops = FALSE
  )

  expect_false(bernoulli$poisson_edges)
  expect_false(bernoulli$allow_self_loops)

  A <- sample_sparse(bernoulli)
  expect_true(all(A@x == 1))
  expect_true(all(Matrix::diag(A) == 0))
})

test_that("directed planted partitions sample graphs", {
  set.seed(27)

  pp <- directed_planted_partition(
    n = 50,
    k = 2,
    within_block = 0.4,
    between_block = 0.05
  )

  A <- sample_sparse(pp)
  expect_equal(dim(A), c(50, 50))

  g <- sample_igraph(pp)
  expect_true(igraph::is_directed(g))
  expect_equal(igraph::vcount(g), 50)
})

test_that("directed_planted_partition() validates its inputs", {
  expect_error(
    directed_planted_partition(n = 10, k = 11, within_block = 0.5, between_block = 0.1),
    "`k` must be less than or equal to `n`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 0, within_block = 0.5, between_block = 0.1),
    "`k` must be a positive integer.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 2.5, k = 2, within_block = 0.5, between_block = 0.1),
    "`n` must be a positive integer.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2),
    "Must specify either `within_block` and `between_block`, or `a` and `b`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, within_block = 0.5),
    "Must specify either `within_block` and `between_block`, or `a` and `b`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, within_block = 0.5, a = 5, between_block = 0.1),
    "Must specify only one of `within_block` and `a`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, within_block = 0.5, between_block = 0.1, b = 1),
    "Must specify only one of `between_block` and `b`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, within_block = c(0.5, 0.4), between_block = 0.1),
    "`within_block` must be a scalar.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, within_block = 1.5, between_block = 0.1),
    "`within_block` (or `a / n`) must be between zero and one.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, within_block = 0.5, between_block = -0.1),
    "`between_block` (or `b / n`) must be between zero and one.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, a = -1, b = 1),
    "`a` must be non-negative.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, a = 1, b = -1),
    "`b` must be non-negative.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(n = 10, k = 2, a = 20, b = 1),
    "`within_block` (or `a / n`) must be between zero and one.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(
      n = 10, k = 2, within_block = 0.5, between_block = 0.1,
      block_sizes = c(5, 5), pi = c(1, 1)
    ),
    "Must specify at most one of `block_sizes` and `pi`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(
      n = 10, k = 2, within_block = 0.5, between_block = 0.1,
      block_sizes = c(5, 3, 2)
    ),
    "Length of `block_sizes` must match `k`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(
      n = 10, k = 2, within_block = 0.5, between_block = 0.1,
      block_sizes = c(5, 6)
    ),
    "Sum of `block_sizes` must equal `n`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(
      n = 10, k = 2, within_block = 0.5, between_block = 0.1,
      pi = c(1, 1, 1)
    ),
    "Length of `pi` must match `k`.",
    fixed = TRUE
  )

  expect_error(
    directed_planted_partition(
      n = 10, k = 2, within_block = 0.5, between_block = 0.1,
      pi = c(1, -1)
    ),
    "All elements of `pi` must be >= 0.",
    fixed = TRUE
  )
})

test_that("print method for directed planted partitions", {
  set.seed(27)

  pp <- directed_planted_partition(
    n = 30,
    k = 3,
    within_block = 0.4,
    between_block = 0.1
  )

  expect_output(print(pp), "Directed Planted Partition Model", fixed = TRUE)
  expect_output(print(pp), "Blocks (k): 3", fixed = TRUE)
  expect_output(print(pp), "Within block edge probability: 0.4", fixed = TRUE)
  expect_output(print(pp), "Between block edge probability: 0.1", fixed = TRUE)
  expect_output(print(pp), "Shared blocks: FALSE", fixed = TRUE)

  shared <- directed_planted_partition(
    n = 30,
    k = 3,
    within_block = 0.4,
    between_block = 0.1,
    shared_blocks = TRUE
  )

  expect_output(print(shared), "Shared blocks: TRUE", fixed = TRUE)
})
