# stored edge probabilities must describe the final model, after rescaling
# by `expected_degree`, `expected_density` and friends

prob_S <- function(model) {
  S <- as.matrix(model$S)
  if (model$poisson_edges) S else 1 - exp(-S)
}

test_that("erdos_renyi() stores the rescaled p", {
  for (poisson_edges in c(TRUE, FALSE)) {
    er <- erdos_renyi(n = 100, expected_degree = 5, poisson_edges = poisson_edges)
    expect_equal(er$p, prob_S(er)[1, 1])
    expect_equal(er$p, 0.05)

    er <- erdos_renyi(n = 100, p = 0.5, expected_density = 0.1, poisson_edges = poisson_edges)
    expect_equal(er$p, prob_S(er)[1, 1])
    expect_equal(er$p, 0.1)

    er <- erdos_renyi(n = 100, p = 0.2, poisson_edges = poisson_edges)
    expect_identical(er$p, 0.2)
  }

  # Poisson rates above one are fine after rescaling
  er <- erdos_renyi(n = 10, expected_degree = 50)
  expect_equal(er$p, 5)
})

test_that("directed_erdos_renyi() stores the rescaled p", {
  for (poisson_edges in c(TRUE, FALSE)) {
    er <- directed_erdos_renyi(n = 100, expected_out_degree = 5, poisson_edges = poisson_edges)
    expect_equal(er$p, prob_S(er)[1, 1])
    expect_equal(er$p, 0.05)

    er <- directed_erdos_renyi(n = 100, p = 0.2, poisson_edges = poisson_edges)
    expect_identical(er$p, 0.2)
  }
})

test_that("erdos_renyi() and directed_erdos_renyi() check p the same way", {
  for (f in list(erdos_renyi, directed_erdos_renyi)) {
    expect_error(f(n = 10, p = 1.5), "strictly between zero and one")
    expect_error(f(n = 10, p = 0), "strictly between zero and one")
    expect_error(f(n = 10, p = -0.1), "strictly between zero and one")
    expect_error(f(n = 10, p = c(0.1, 0.2)), "strictly between zero and one")
  }
})

test_that("planted_partition() stores the rescaled edge probabilities", {
  for (poisson_edges in c(TRUE, FALSE)) {
    pp <- planted_partition(
      n = 100, k = 2, within_block = .5, between_block = .1,
      expected_degree = 5, poisson_edges = poisson_edges
    )
    S <- prob_S(pp)
    expect_equal(pp$within_block, S[1, 1])
    expect_equal(pp$between_block, S[1, 2])
    # the ratio is preserved
    expect_equal(pp$within_block / pp$between_block, 5)
  }

  # without rescaling the given values are kept exactly
  pp <- planted_partition(n = 100, k = 2, within_block = .5, between_block = .1)
  expect_identical(c(pp$within_block, pp$between_block), c(.5, .1))

  # Poisson rates above one are fine after rescaling
  pp <- planted_partition(
    n = 10, k = 2, within_block = .5, between_block = .1, expected_degree = 30
  )
  expect_gt(pp$within_block, 1)
  expect_equal(pp$within_block, as.matrix(pp$S)[1, 1])
})

test_that("directed_planted_partition() stores the rescaled edge probabilities", {
  for (poisson_edges in c(TRUE, FALSE)) {
    pp <- directed_planted_partition(
      n = 100, k = 2, within_block = .5, between_block = .1,
      expected_out_degree = 5, poisson_edges = poisson_edges
    )
    S <- prob_S(pp)
    expect_equal(pp$within_block, S[1, 1])
    expect_equal(pp$between_block, S[1, 2])
    expect_equal(pp$within_block / pp$between_block, 5)
  }
})
