test_that("Bernoulli models with an edge probability of 1 error clearly", {
  msg <- "must be less than 1 for Bernoulli graphs"

  expect_error(
    planted_partition(
      n = 10, k = 2, within_block = 1, between_block = 0.1,
      poisson_edges = FALSE
    ),
    msg
  )
  expect_error(
    sbm(n = 10, B = matrix(c(1, .1, .1, .5), 2), poisson_edges = FALSE),
    msg
  )
  expect_error(
    undirected_factor_model(matrix(1, 5, 1), matrix(1), poisson_edges = FALSE),
    msg
  )
  expect_error(
    directed_factor_model(
      matrix(1, 5, 1), matrix(1), matrix(1, 4, 1),
      poisson_edges = FALSE
    ),
    msg
  )
  expect_error(
    directed_planted_partition(
      n = 10, k = 2, within_block = 1, between_block = 0.1,
      poisson_edges = FALSE
    ),
    msg
  )
  # rescaling to a probability of 1
  expect_error(
    erdos_renyi(n = 10, expected_degree = 10, poisson_edges = FALSE),
    msg
  )

  # Poisson rates of 1 are fine
  m <- sbm(n = 10, B = matrix(c(1, .1, .1, .5), 2))
  expect_true(all(is.finite(as.matrix(m$S))))
})

test_that("poisson_edges and allow_self_loops must be TRUE or FALSE", {
  for (bad in list("yes", NA, c(TRUE, FALSE), 1)) {
    expect_error(
      undirected_factor_model(matrix(1, 5, 1), matrix(.1), poisson_edges = bad),
      "`poisson_edges` must be `TRUE` or `FALSE`.",
      fixed = TRUE
    )
    expect_error(
      directed_factor_model(
        matrix(1, 5, 1), matrix(.1), matrix(1, 4, 1),
        poisson_edges = bad
      ),
      "`poisson_edges` must be `TRUE` or `FALSE`.",
      fixed = TRUE
    )
    expect_error(
      sbm(n = 10, B = diag(2) * .5, allow_self_loops = bad),
      "`allow_self_loops` must be `TRUE` or `FALSE`.",
      fixed = TRUE
    )
  }
})
