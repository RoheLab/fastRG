test_that("directed_dcsbm() errors when k_out / k_in contradict B", {
  B <- matrix(0.1, 2, 3)

  expect_error(
    directed_dcsbm(n = 10, B = B, k_out = 5, k_in = 3),
    "`k_out` must equal `nrow(B)`",
    fixed = TRUE
  )
  expect_error(
    directed_dcsbm(n = 10, B = B, k_out = 2, k_in = 5),
    "`k_in` must equal `ncol(B)`",
    fixed = TRUE
  )

  # agreeing values are still fine, as documented
  m <- suppressMessages(directed_dcsbm(n = 10, B = B, k_out = 2, k_in = 3))
  expect_equal(c(m$k1, m$k2), c(2, 3))
})

test_that("mmsbm() rejects redundant arguments and degenerate k", {
  B <- diag(2) * 0.5 + 0.1

  expect_error(
    mmsbm(theta = rep(1, 10), k = 3, B = B),
    "Must specify only one of `B` and `k`, not both.",
    fixed = TRUE
  )
  expect_error(
    mmsbm(n = 100, theta = rep(1, 10), B = B),
    "Must specify only one of `n` and `theta`, not both.",
    fixed = TRUE
  )
  expect_error(
    mmsbm(theta = rep(1, 10), B = matrix(0.5)),
    "at least two blocks"
  )
  expect_error(
    mmsbm(theta = rep(1, 2), B = diag(3) + 0.1),
    "`force_pure = TRUE` needs at least as many nodes as blocks",
    fixed = TRUE
  )

  # without force_pure, n < k is fine
  m <- mmsbm(theta = rep(1, 2), B = diag(3) + 0.1, force_pure = FALSE)
  expect_equal(m$n, 2)
})

test_that("overlapping_sbm() errors clearly when force_pure and n < k", {
  expect_error(
    overlapping_sbm(n = 2, B = diag(3) * 0.5 + 0.1),
    "`force_pure = TRUE` needs at least as many nodes as blocks",
    fixed = TRUE
  )
})

test_that("chung_lu() errors when both n and theta are given", {
  expect_error(
    chung_lu(n = 100, theta = rep(1, 10)),
    "Must specify only one of `n` and `theta`, not both.",
    fixed = TRUE
  )
})
