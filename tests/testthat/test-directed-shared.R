# shared blocks and shared degree heterogeneity in directed_dcsbm() (#40)

block_of_row <- function(X) {
  X <- methods::as(X, "TsparseMatrix")
  out <- integer(nrow(X))
  out[X@i + 1] <- X@j + 1
  out
}

test_that("explicit shared_* = FALSE matches the default", {
  set.seed(27)
  default <- suppressMessages(
    directed_dcsbm(n = 30, k_out = 3, k_in = 2, expected_density = 0.1)
  )

  set.seed(27)
  explicit <- suppressMessages(
    directed_dcsbm(
      n = 30,
      k_out = 3,
      k_in = 2,
      expected_density = 0.1,
      shared_blocks = FALSE,
      shared_theta = FALSE
    )
  )

  expect_identical(default, explicit)
  expect_false(default$shared_blocks)
  expect_false(default$shared_theta)
})

test_that("shared_blocks uses one membership vector for X and Y", {
  B <- matrix(1:9, nrow = 3, ncol = 3)
  theta_out <- stats::rexp(60) + 1
  theta_in <- stats::rexp(60) + 1

  for (sort_nodes in c(TRUE, FALSE)) {
    set.seed(27)
    model <- directed_dcsbm(
      theta_out = theta_out,
      theta_in = theta_in,
      B = B,
      shared_blocks = TRUE,
      sort_nodes = sort_nodes
    )

    expect_true(model$shared_blocks)
    expect_identical(model$z_out, model$z_in)
    expect_identical(model$pi_out, model$pi_in)

    # row i of X and row i of Y are the same node, in the same block
    expect_equal(block_of_row(model$X), as.integer(model$z_out))
    expect_equal(block_of_row(model$Y), as.integer(model$z_in))

    # stored theta matches the rows of X and Y
    expect_equal(Matrix::rowSums(model$X), model$theta_out, ignore_attr = TRUE)
    expect_equal(Matrix::rowSums(model$Y), model$theta_in, ignore_attr = TRUE)

    if (sort_nodes) {
      expect_false(is.unsorted(model$z_out))
    }
  }
})

test_that("shared_blocks keeps node identity under sorting", {
  theta_out <- as.numeric(1:20)
  theta_in <- as.numeric(20:1)

  set.seed(27)
  model <- directed_dcsbm(
    theta_out = theta_out,
    theta_in = theta_in,
    B = diag(2),
    shared_blocks = TRUE,
    force_identifiability = FALSE
  )

  # theta_out + theta_in is 21 for every node, which is only preserved if
  # rows of X and rows of Y are permuted together
  expect_equal(model$theta_out + model$theta_in, rep(21, 20))
  expect_equal(
    Matrix::rowSums(model$X) + Matrix::rowSums(model$Y),
    rep(21, 20),
    ignore_attr = TRUE
  )
})

test_that("shared_blocks permutes B and block sizes consistently", {
  B <- matrix(1:9, nrow = 3, ncol = 3)

  set.seed(27)
  model <- directed_dcsbm(
    theta_out = rep(1, 100),
    theta_in = rep(1, 100),
    B = B,
    block_sizes_in = c(50, 20, 30),
    shared_blocks = TRUE,
    force_identifiability = FALSE
  )

  expect_identical(model$z_out, model$z_in)
  expect_equal(as.vector(table(model$z_out)), c(20, 30, 50))
  expect_equal(unname(as.matrix(model$S)), B[c(2, 3, 1), c(2, 3, 1)])

  set.seed(27)
  unsorted <- directed_dcsbm(
    theta_out = rep(1, 100),
    theta_in = rep(1, 100),
    B = B,
    block_sizes_out = c(50, 20, 30),
    block_sizes_in = c(50, 20, 30),
    shared_blocks = TRUE,
    sort_nodes = FALSE
  )

  expect_identical(unsorted$z_out, unsorted$z_in)
  expect_equal(as.vector(table(unsorted$z_out)), c(50, 20, 30))
  expect_equal(unname(as.matrix(unsorted$S)), B)
})

test_that("shared_blocks accepts pi from either side, or the default", {
  B <- matrix(0.5, nrow = 3, ncol = 3)
  theta <- rep(1, 30)

  set.seed(27)
  from_out <- directed_dcsbm(
    theta_out = theta,
    theta_in = theta,
    B = B,
    pi_out = c(1, 2, 3),
    shared_blocks = TRUE
  )

  set.seed(27)
  from_in <- directed_dcsbm(
    theta_out = theta,
    theta_in = theta,
    B = B,
    pi_in = c(1, 2, 3),
    shared_blocks = TRUE
  )

  set.seed(27)
  from_both <- directed_dcsbm(
    theta_out = theta,
    theta_in = theta,
    B = B,
    pi_out = c(1, 2, 3),
    pi_in = c(1, 2, 3) / 6,
    shared_blocks = TRUE
  )

  expect_identical(from_out, from_in)
  expect_equal(from_out, from_both)
  expect_equal(from_out$pi_out, c(1, 2, 3) / 6)

  set.seed(27)
  default_pi <- directed_dcsbm(
    theta_out = theta,
    theta_in = theta,
    B = B,
    shared_blocks = TRUE
  )

  expect_equal(default_pi$pi_out, rep(1 / 3, 3))
  expect_identical(default_pi$z_out, default_pi$z_in)
})

test_that("shared_blocks needs only one of k_out and k_in", {
  set.seed(27)
  model <- suppressMessages(
    directed_dcsbm(n = 20, k_in = 4, shared_blocks = TRUE)
  )

  expect_equal(dim(model$S), c(4, 4))
  expect_identical(model$z_out, model$z_in)
})

test_that("shared_blocks validates its inputs", {
  theta <- rep(1, 10)

  expect_error(
    directed_dcsbm(
      theta_out = theta,
      theta_in = theta,
      B = matrix(0.5, 2, 3),
      shared_blocks = TRUE
    ),
    "`B` must be a square matrix when `shared_blocks = TRUE`.",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(
      theta_out = theta,
      theta_in = theta,
      k_out = 2,
      k_in = 3,
      shared_blocks = TRUE
    ),
    "`k_out` and `k_in` must be equal when `shared_blocks = TRUE`.",
    fixed = TRUE
  )

  B <- matrix(0.5, 2, 2)

  expect_error(
    directed_dcsbm(
      theta_out = theta,
      theta_in = theta,
      B = B,
      pi_out = c(1, 2),
      pi_in = c(2, 1),
      shared_blocks = TRUE
    ),
    "`pi_out` and `pi_in` must agree",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(
      theta_out = theta,
      theta_in = theta,
      B = B,
      block_sizes_out = c(4, 6),
      block_sizes_in = c(6, 4),
      shared_blocks = TRUE
    ),
    "`block_sizes_out` and `block_sizes_in` must be identical",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(
      theta_out = theta,
      theta_in = theta,
      B = B,
      block_sizes_out = c(4, 6),
      pi_in = c(1, 1),
      shared_blocks = TRUE
    ),
    "must specify at most one of block sizes",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(
      theta_out = theta,
      theta_in = theta,
      B = B,
      shared_blocks = NA
    ),
    "`shared_blocks` must be `TRUE` or `FALSE`.",
    fixed = TRUE
  )
})

test_that("shared_theta uses one degree heterogeneity vector", {
  B <- matrix(1:6, nrow = 2, ncol = 3)
  theta <- stats::rexp(40) + 1

  for (sort_nodes in c(TRUE, FALSE)) {
    set.seed(27)
    model <- directed_dcsbm(
      theta_in = theta,
      B = B,
      shared_theta = TRUE,
      force_identifiability = FALSE,
      sort_nodes = sort_nodes
    )

    expect_true(model$shared_theta)
    expect_false(model$shared_blocks)
    expect_identical(model$theta_out, model$theta_in)
    expect_setequal(model$theta_out, theta)

    # row i of X and row i of Y belong to the same node
    expect_equal(Matrix::rowSums(model$X), model$theta_out, ignore_attr = TRUE)
    expect_equal(Matrix::rowSums(model$Y), model$theta_in, ignore_attr = TRUE)
    expect_equal(block_of_row(model$X), as.integer(model$z_out))
    expect_equal(block_of_row(model$Y), as.integer(model$z_in))

    if (sort_nodes) {
      expect_false(is.unsorted(model$z_out))
    } else {
      expect_identical(model$theta_out, theta)
    }
  }

  set.seed(27)
  from_out <- directed_dcsbm(
    theta_out = theta,
    B = B,
    shared_theta = TRUE
  )

  set.seed(27)
  from_both <- directed_dcsbm(
    theta_out = theta,
    theta_in = theta,
    B = B,
    shared_theta = TRUE
  )

  expect_identical(from_out, from_both)
})

test_that("shared_theta generates a single random vector from n", {
  set.seed(27)
  expect_message(
    model <- directed_dcsbm(
      n = 25,
      B = diag(2),
      shared_theta = TRUE,
      force_identifiability = FALSE,
      sort_nodes = FALSE
    ),
    "LogNormal(2, 1)",
    fixed = TRUE
  )

  set.seed(27)
  expected_theta <- stats::rlnorm(25, meanlog = 2, sdlog = 1)

  expect_identical(model$theta_out, expected_theta)
  expect_identical(model$theta_in, expected_theta)
})

test_that("shared_theta with force_identifiability normalizes within blocks", {
  set.seed(27)
  model <- directed_dcsbm(
    theta_out = stats::rexp(50) + 1,
    B = matrix(1, 2, 3),
    shared_theta = TRUE,
    force_identifiability = TRUE
  )

  # normalized separately within outgoing / incoming blocks, so theta_out
  # and theta_in are proportional within each pair of blocks
  ratio <- model$theta_out / model$theta_in
  cell <- interaction(model$z_out, model$z_in, drop = TRUE)
  ratio_spread <- tapply(ratio, cell, function(r) diff(range(r)))
  expect_true(all(ratio_spread < 1e-12))

  # and equal when blocks are shared too
  set.seed(27)
  both <- directed_dcsbm(
    theta_out = stats::rexp(50) + 1,
    B = diag(3),
    shared_theta = TRUE,
    shared_blocks = TRUE
  )

  expect_identical(both$theta_out, both$theta_in)
  expect_identical(both$z_out, both$z_in)
  expect_equal(
    as.vector(tapply(both$theta_out, both$z_out, sum)),
    rep(1, 3)
  )
})

test_that("shared_theta validates its inputs", {
  B <- diag(2)

  expect_error(
    directed_dcsbm(
      theta_out = c(1, 2, 3),
      theta_in = c(1, 2, 4),
      B = B,
      shared_theta = TRUE
    ),
    "`theta_out` and `theta_in` must be identical",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(
      theta_out = c(1, 2, 3),
      theta_in = c(1, 2),
      B = B,
      shared_theta = TRUE
    ),
    "`theta_out` and `theta_in` must be identical",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(n = 3, theta_out = c(1, 2, 3), B = B, shared_theta = TRUE),
    "Must specify only one of `n` and `theta_out` / `theta_in`, not both.",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(B = B, shared_theta = TRUE),
    "Must specify either `n`, or `theta_out` or `theta_in`.",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(n = 3, theta_out = c(1, 2, 3), B = B),
    "Must specify either `n`, or both `theta_out` and `theta_in` together",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(theta_out = c(1, 2), B = B),
    "Must specify either `n` or both `theta_out` and `theta_in` together.",
    fixed = TRUE
  )

  expect_error(
    directed_dcsbm(n = 3, B = B, shared_theta = "yes"),
    "`shared_theta` must be `TRUE` or `FALSE`.",
    fixed = TRUE
  )
})

test_that("print method shows shared structure", {
  set.seed(27)

  shared <- directed_dcsbm(
    theta_out = rep(1, 10),
    B = diag(2),
    shared_blocks = TRUE,
    shared_theta = TRUE
  )

  expect_output(print(shared), "Shared blocks: TRUE", fixed = TRUE)
  expect_output(
    print(shared),
    "Shared degree heterogeneity: TRUE",
    fixed = TRUE
  )
  expect_output(print(shared), "Outgoing Blocks (k_out): 2", fixed = TRUE)

  separate <- directed_dcsbm(
    theta_out = rep(1, 10),
    theta_in = rep(1, 10),
    B = matrix(1, 2, 3)
  )

  expect_output(print(separate), "Shared blocks: FALSE", fixed = TRUE)
  expect_output(
    print(separate),
    "Shared degree heterogeneity: FALSE",
    fixed = TRUE
  )
  expect_output(print(separate), "Incoming Blocks (k_in): 3", fixed = TRUE)
})
