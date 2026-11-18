test_that("dcsbm() errors when both `pi` and `block_sizes` are given", {
  expect_error(
    dcsbm(theta = rep(1, 10), B = diag(2), block_sizes = c(5, 5), pi = c(.9, .1)),
    "only one of `block_sizes` and `pi`"
  )
  expect_error(
    sbm(n = 10, B = diag(2), block_sizes = c(5, 5), pi = c(.9, .1)),
    "only one of `block_sizes` and `pi`"
  )
})

test_that("dcsbm() validates `block_sizes`", {
  # non-integer sizes used to silently drop nodes
  expect_error(
    dcsbm(theta = 1:5 + 0, B = diag(2), block_sizes = c(2.5, 2.5)),
    "non-negative integers"
  )
  expect_error(
    dcsbm(theta = rep(1, 10), B = diag(2), block_sizes = c(3, 3, 4)),
    "Length of `block_sizes` must match"
  )
  expect_error(
    dcsbm(theta = rep(1, 10), B = diag(2), block_sizes = c(-1, 11)),
    "non-negative integers"
  )
  expect_error(
    dcsbm(theta = rep(1, 10), B = diag(2), block_sizes = c(4, 4)),
    "Sum of `block_sizes` must equal"
  )
})

test_that("dcsbm() uses `pi` or `block_sizes` when given alone", {
  set.seed(1)
  m <- dcsbm(theta = rep(1, 1000), B = diag(2), pi = c(.9, .1), sort_nodes = FALSE)
  # exact binomial bound on the size of the first block
  expect_lt(abs(sum(m$z == "block1") - 900) / sqrt(1000 * .9 * .1), 5)

  m <- dcsbm(theta = rep(1, 10), B = diag(2), block_sizes = c(3, 7), sort_nodes = FALSE)
  expect_equal(as.vector(table(m$z)), c(3, 7))
  expect_equal(m$pi, c(.3, .7))

  m <- dcsbm(theta = rep(1, 11), B = diag(3))
  expect_equal(sort(as.vector(table(m$z))), c(3, 4, 4))
})
