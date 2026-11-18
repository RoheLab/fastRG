test_that("planted_partition() has its own subclass", {
  set.seed(27)

  pp <- planted_partition(
    n = 60,
    k = 3,
    within_block = 0.3,
    between_block = 0.05
  )

  expect_s3_class(
    pp,
    c(
      "undirected_planted_partition",
      "undirected_sbm",
      "undirected_dcsbm",
      "undirected_factor_model"
    ),
    exact = TRUE
  )

  expect_equal(pp$within_block, 0.3)
  expect_equal(pp$between_block, 0.05)

  expected_B <- matrix(0.05, 3, 3)
  diag(expected_B) <- 0.3
  expect_equal(unname(as.matrix(pp$S)), expected_B)
})

test_that("planted_partition() samples exactly like the equivalent sbm()", {
  B <- matrix(0.05, 4, 4)
  diag(B) <- 0.3

  set.seed(27)
  pp <- planted_partition(
    n = 80,
    k = 4,
    within_block = 0.3,
    between_block = 0.05,
    expected_degree = 5
  )

  set.seed(27)
  model <- sbm(n = 80, B = B, expected_degree = 5)

  expect_identical(pp$z, model$z)
  expect_identical(pp$X, model$X)
  expect_identical(pp$S, model$S)
})

test_that("planted partitions print via their own methods", {
  set.seed(27)

  undirected <- planted_partition(
    n = 30,
    k = 3,
    within_block = 0.4,
    between_block = 0.1
  )

  undirected_output <- capture.output(print(undirected))
  expect_identical(
    undirected_output[1],
    "Undirected Planted Partition Model"
  )
  expect_true(
    "Within block edge probability: 0.4 " %in% undirected_output
  )
  expect_true(
    "Between block edge probability: 0.1 " %in% undirected_output
  )

  directed <- directed_planted_partition(
    n = 30,
    k = 3,
    within_block = 0.4,
    between_block = 0.1
  )

  directed_output <- capture.output(print(directed))
  expect_identical(directed_output[1], "Directed Planted Partition Model")
  expect_true(
    "Within block edge probability: 0.4 " %in% directed_output
  )
})

test_that("planted partition validators check B", {
  set.seed(27)

  pp <- planted_partition(
    n = 30,
    k = 3,
    within_block = 0.4,
    between_block = 0.1
  )

  pp$S[1, 2] <- 0.2
  expect_error(
    validate_undirected_planted_partition(pp),
    "The off-diagonal elements of `B` must be constant.",
    fixed = TRUE
  )

  dpp <- directed_planted_partition(
    n = 30,
    k = 3,
    within_block = 0.4,
    between_block = 0.1
  )

  dpp$S[2, 2] <- 0.9
  expect_error(
    validate_directed_planted_partition(dpp),
    "The diagonal of `B` must be constant.",
    fixed = TRUE
  )

  expect_error(
    planted_partition(n = 30, k = 3, within_block = 1.5, between_block = 0.1),
    "`within_block` must be between zero and one.",
    fixed = TRUE
  )
})

test_that("planted_partition() requires exactly one way to give each probability", {
  expect_error(
    planted_partition(n = 10, k = 2, within_block = 0.5, a = 2, between_block = 0.1),
    "only one of `within_block` and `a`"
  )
  expect_error(
    planted_partition(n = 10, k = 2, within_block = 0.5, between_block = 0.1, b = 2),
    "only one of `between_block` and `b`"
  )
  expect_error(
    planted_partition(n = 10, k = 2),
    "Must specify either `within_block` and `between_block`, or `a` and `b`"
  )
  expect_error(
    planted_partition(n = 10, k = 2, within_block = 0.5),
    "Must specify either `within_block` and `between_block`, or `a` and `b`"
  )
  expect_error(
    planted_partition(n = 10, k = 2, within_block = c(0.5, 0.2), between_block = 0.1),
    "`within_block` must be a scalar"
  )

  # mixing the two parameterizations is fine
  pp <- planted_partition(n = 10, k = 2, within_block = 0.5, b = 2)
  expect_equal(pp$between_block, 0.2)
})
