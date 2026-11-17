sparse_triplets <- function(A) {
  s <- Matrix::summary(methods::as(A, "TsparseMatrix"))
  list(i = as.integer(s$i), j = as.integer(s$j), x = as.numeric(s$x))
}

sample_reference_case <- function(...) {
  set.seed(27)
  model <- suppressMessages(directed_dcsbm(...))
  list(model = model, A = sample_sparse(model))
}

expect_reference_case <- function(case, z_out, z_in, theta_out, theta_in, S, i, j, x) {
  model <- case$model
  expect_equal(as.integer(model$z_out), z_out)
  expect_equal(as.integer(model$z_in), z_in)
  expect_equal(unname(model$theta_out), theta_out, tolerance = 1e-6)
  expect_equal(unname(model$theta_in), theta_in, tolerance = 1e-6)
  expect_equal(unname(as.matrix(model$S)), S, tolerance = 1e-6)
  expect_equal(sparse_triplets(case$A), list(i = i, j = j, x = x))
}

# reference values pin seeded output so that accidental changes to sampling
# behavior (or RNG order) are caught.
#
# - the default-case pins were regenerated when the directed default changed
#   from sampling memberships with balanced `pi` to exactly balanced block
#   sizes (matching `dcsbm()`), which intentionally changes seeded output
#   of the defaults. they also reflect the theta-order fix (theta in node
#   order).
# - the explicit `pi` pin was generated with fastRG 0.4.0.9000 (commit
#   30679ae), before `block_sizes_out` / `block_sizes_in` were added, and is
#   unchanged

test_that("directed_dcsbm() default sampling is pinned (sorted)", {
  case <- sample_reference_case(
    n = 15,
    k_out = 3,
    k_in = 4,
    expected_density = 0.2
  )

  expect_reference_case(
    case,
    z_out = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
    z_in = c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L),
    theta_out = c(
      0.617259, 0.2880123, 0.04267424, 0.03071212, 0.02134236, 0.5541946,
      0.208573, 0.08806888, 0.08315635, 0.06600721, 0.619145, 0.1508151,
      0.1333359, 0.06108368, 0.03562034
    ),
    theta_in = c(
      0.3888567, 0.3887706, 0.2223727, 0.4261955, 0.3017147, 0.1742412,
      0.09784863, 0.4608853, 0.3938223, 0.1242239, 0.02106845, 0.5364644,
      0.2054845, 0.2038565, 0.05419462
    ),
    S = structure(
      c(
        4.604851, 6.235598, 1.309489, 3.287976, 3.37821,
        1.262552, 4.444437, 7.243754, 6.082634, 0.4312557, 4.136733,
        2.58251
      ),
      dim = 3:4
    ),
    i = c(
      1L, 2L, 6L, 7L, 8L, 13L, 1L, 3L, 6L, 9L, 10L, 1L, 2L, 5L, 6L,
      12L, 1L, 6L, 1L, 2L, 6L, 7L, 8L, 11L, 6L, 9L, 2L, 1L, 6L, 9L,
      10L, 11L, 12L, 13L, 6L, 7L, 10L, 13L, 6L, 7L, 11L, 7L, 11L, 6L,
      11L
    ),
    j = c(
      1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L,
      3L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 7L, 8L, 8L, 8L, 8L,
      8L, 8L, 8L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 12L, 12L, 13L, 14L
    ),
    x = c(
      1, 2, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
      3, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 3,
      1, 1, 1, 1
    )
  )
})

test_that("directed_dcsbm() default sampling is pinned (unsorted)", {
  case <- sample_reference_case(
    n = 15,
    k_out = 3,
    k_in = 4,
    expected_density = 0.2,
    sort_nodes = FALSE
  )

  expect_reference_case(
    case,
    z_out = c(3L, 1L, 3L, 1L, 1L, 3L, 2L, 3L, 2L, 3L, 2L, 1L, 2L, 1L, 2L),
    z_in = c(3L, 3L, 2L, 3L, 3L, 1L, 4L, 1L, 4L, 1L, 4L, 2L, 1L, 2L, 2L),
    theta_out = c(
      0.5182694, 0.5310611, 0.03583059, 0.03935283, 0.05662958, 0.1033968,
      0.06703823, 0.2448739, 0.562851, 0.09762926, 0.01841655, 0.1750073,
      0.3201125, 0.1979492, 0.03158169
    ),
    theta_in = c(
      0.1187938, 0.06794876, 0.05971175, 0.3370907, 0.4761668, 0.1813899,
      0.2049101, 0.3516366, 0.6044328, 0.01607437, 0.1906571, 0.06413791,
      0.4508992, 0.6348915, 0.2412588
    ),
    S = structure(
      c(
        3.287976, 3.37821, 1.262552, 4.444437, 7.243754,
        6.082634, 0.4312557, 4.136733, 2.58251, 4.604851, 6.235598, 1.309489
      ),
      dim = 3:4
    ),
    i = c(
      15L, 1L, 9L, 6L, 1L, 9L, 11L, 13L, 12L, 1L, 9L, 2L, 4L, 13L,
      1L, 9L, 13L, 2L, 7L, 13L, 2L, 12L, 2L, 8L, 9L, 10L, 12L, 13L,
      14L, 1L, 2L, 5L, 7L, 8L, 9L, 13L, 14L, 1L, 2L, 6L, 9L, 13L, 14L,
      15L
    ),
    j = c(
      1L, 3L, 3L, 4L, 5L, 5L, 5L, 5L, 6L, 7L, 7L, 8L, 8L, 8L, 9L,
      9L, 9L, 11L, 11L, 11L, 12L, 12L, 13L, 13L, 13L, 13L, 13L, 13L,
      13L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 15L, 15L, 15L, 15L,
      15L, 15L, 15L
    ),
    x = c(
      1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 6, 2, 1, 1, 1, 1, 1,
      1, 1, 1
    )
  )
})

test_that("directed_dcsbm() defaults to balanced block sizes", {
  set.seed(27)

  model <- suppressMessages(
    directed_dcsbm(n = 23, k_out = 4, k_in = 5, sort_nodes = FALSE)
  )

  # 23 = 6 + 6 + 6 + 5 and 23 = 5 + 5 + 5 + 4 + 4
  expect_equal(as.vector(table(model$z_out)), c(6, 6, 6, 5))
  expect_equal(as.vector(table(model$z_in)), c(5, 5, 5, 4, 4))
  expect_equal(model$pi_out, c(6, 6, 6, 5) / 23)
  expect_equal(model$pi_in, c(5, 5, 5, 4, 4) / 23)

  set.seed(27)

  sorted <- suppressMessages(
    directed_dcsbm(n = 23, k_out = 4, k_in = 5)
  )

  # sorting relabels blocks in increasing order of size
  expect_equal(as.vector(table(sorted$z_out)), c(5, 6, 6, 6))
  expect_equal(as.vector(table(sorted$z_in)), c(4, 4, 5, 5, 5))
})

test_that("directed_dcsbm() sampling with explicit pi is unchanged", {
  case <- sample_reference_case(
    theta_out = rep(1:3, 5),
    theta_in = rep(2, 15),
    B = matrix(1:6, 2, 3),
    pi_out = c(3, 1),
    pi_in = c(1, 2, 3),
    expected_density = 0.2
  )

  expect_reference_case(
    case,
    z_out = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L),
    z_in = c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L),
    theta_out = c(
      0.4285714, 0.2857143, 0.1428571, 0.1428571, 0.1304348, 0.1304348,
      0.1304348, 0.1304348, 0.08695652, 0.08695652, 0.08695652, 0.08695652,
      0.04347826, 0.04347826, 0.04347826
    ),
    theta_in = c(
      0.3333333, 0.3333333, 0.3333333, 0.25, 0.25, 0.25, 0.25, 0.125,
      0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125
    ),
    S = structure(
      c(4.285714, 2.142857, 8.571429, 6.428571, 12.85714, 10.71429),
      dim = 2:3
    ),
    i = c(
      1L, 2L, 1L, 3L, 7L, 8L, 12L, 1L, 7L, 2L, 4L, 5L, 15L, 1L, 7L,
      12L, 1L, 2L, 5L, 10L, 6L, 8L, 2L, 4L, 2L, 9L, 2L, 11L, 3L, 1L,
      3L
    ),
    j = c(
      1L, 1L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 6L, 7L, 7L, 7L, 8L, 8L,
      8L, 9L, 9L, 9L, 9L, 10L, 10L, 11L, 11L, 12L, 12L, 13L, 13L, 14L,
      15L, 15L
    ),
    x = c(
      1, 1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
      1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1
    )
  )
})

test_that("directed_dcsbm() respects block_sizes_out and block_sizes_in", {
  set.seed(27)

  B <- matrix(1:12, nrow = 3, ncol = 4)

  for (sort_nodes in c(TRUE, FALSE)) {
    model <- directed_dcsbm(
      theta_out = rep(1, 100),
      theta_in = rep(1, 100),
      B = B,
      block_sizes_out = c(50, 20, 30),
      block_sizes_in = c(10, 40, 25, 25),
      sort_nodes = sort_nodes
    )

    # block sizes are exact, whatever the ordering of the blocks
    expect_equal(sort(as.vector(table(model$z_out))), c(20, 30, 50))
    expect_equal(sort(as.vector(table(model$z_in))), c(10, 25, 25, 40))

    # stored block probabilities match observed block sizes
    expect_equal(as.vector(table(model$z_out)) / 100, model$pi_out)
    expect_equal(as.vector(table(model$z_in)) / 100, model$pi_in)

    # X and Y have the right number of nodes in each block
    expect_equal(
      unname(Matrix::colSums(model$X != 0)),
      as.vector(table(model$z_out))
    )
    expect_equal(
      unname(Matrix::colSums(model$Y != 0)),
      as.vector(table(model$z_in))
    )
  }
})

test_that("sort_nodes reorders B consistently with block_sizes", {
  set.seed(27)

  B <- matrix(1:12, nrow = 3, ncol = 4)

  model <- directed_dcsbm(
    theta_out = rep(1, 100),
    theta_in = rep(1, 100),
    B = B,
    block_sizes_out = c(50, 20, 30),
    block_sizes_in = c(10, 40, 35, 15),
    sort_nodes = TRUE,
    force_identifiability = FALSE
  )

  expect_true(model$sorted)
  expect_false(is.unsorted(model$z_out))
  expect_false(is.unsorted(model$z_in))

  # blocks are relabeled in increasing order of size
  expect_equal(as.vector(table(model$z_out)), c(20, 30, 50))
  expect_equal(as.vector(table(model$z_in)), c(10, 15, 35, 40))

  # and the rows / columns of B travel with their blocks
  expect_equal(unname(as.matrix(model$S)), B[c(2, 3, 1), c(1, 4, 3, 2)])

  unsorted <- directed_dcsbm(
    theta_out = rep(1, 100),
    theta_in = rep(1, 100),
    B = B,
    block_sizes_out = c(50, 20, 30),
    block_sizes_in = c(10, 40, 35, 15),
    sort_nodes = FALSE,
    force_identifiability = FALSE
  )

  expect_equal(as.vector(table(unsorted$z_out)), c(50, 20, 30))
  expect_equal(as.vector(table(unsorted$z_in)), c(10, 40, 35, 15))
  expect_equal(unname(as.matrix(unsorted$S)), B)
})

test_that("block sizes can be mixed with pi and with random B", {
  set.seed(27)

  model <- suppressMessages(
    directed_dcsbm(
      n = 60,
      k_out = 2,
      k_in = 3,
      block_sizes_out = c(45, 15),
      pi_in = c(1, 1, 1)
    )
  )

  expect_equal(as.vector(table(model$z_out)), c(15, 45))
  expect_equal(nlevels(model$z_in), 3)
})

test_that("block sizes may contain empty blocks", {
  set.seed(27)

  model <- directed_dcsbm(
    theta_out = rep(1, 1),
    theta_in = rep(1, 1),
    B = matrix(1, 2, 2),
    block_sizes_out = c(0, 1),
    block_sizes_in = c(1, 0),
    sort_nodes = FALSE
  )

  expect_equal(length(model$z_out), 1)
  expect_equal(as.integer(model$z_out), 2L)
  expect_equal(as.integer(model$z_in), 1L)
})

test_that("directed_dcsbm() validates block sizes", {
  B <- matrix(0.5, nrow = 2, ncol = 3)
  theta <- rep(1, 10)

  dcsbm_with <- function(...) {
    directed_dcsbm(theta_out = theta, theta_in = theta, B = B, ...)
  }

  expect_error(
    dcsbm_with(block_sizes_out = c(5, 5), pi_out = c(1, 1)),
    "Must specify at most one of `block_sizes_out` and `pi_out`.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_in = c(5, 3, 2), pi_in = c(1, 1, 1)),
    "Must specify at most one of `block_sizes_in` and `pi_in`.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_out = c(5, 3, 2)),
    "Length of `block_sizes_out` must match number of rows in `B`.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_in = c(5, 5)),
    "Length of `block_sizes_in` must match number of columns in `B`.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_out = c(5, 6)),
    "Sum of `block_sizes_out` must equal `n` or `length(theta_out)`.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_in = c(5, 3, 1)),
    "Sum of `block_sizes_in` must equal `n` or `length(theta_in)`.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_out = c(11, -1)),
    "All elements of `block_sizes_out` must be non-negative integers.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_in = c(4.5, 4.5, 1)),
    "All elements of `block_sizes_in` must be non-negative integers.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_out = c(NA, 10)),
    "All elements of `block_sizes_out` must be non-negative integers.",
    fixed = TRUE
  )

  expect_error(
    dcsbm_with(block_sizes_out = c("5", "5")),
    "All elements of `block_sizes_out` must be non-negative integers.",
    fixed = TRUE
  )
})
