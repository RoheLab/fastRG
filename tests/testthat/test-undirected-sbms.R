test_that("SBMs don't drop isolated nodes", {
  set.seed(27)

  n <- 10
  k <- 10
  pi <- rep(1, k) / k

  B <- diag(rep(0.5, k))

  sbm <- sbm(n = n, pi = pi, B = B)

  A <- sample_sparse(sbm)

  expect_equal(ncol(A), 10)
  expect_equal(nrow(A), 10)
})

test_that("rank 1 sbms sample", {
  set.seed(27)

  n <- 10
  k <- 1
  pi <- rep(1, k) / k

  B <- matrix(0.5)

  sbm <- sbm(n = n, pi = pi, B = B)

  expect_silent(A <- sample_sparse(sbm))
})

test_that("validate_undirected_sbm() rejects theta that varies within a block", {
  set.seed(1)
  m <- sbm(n = 10, B = diag(2) * 0.5 + 0.1)
  expect_identical(validate_undirected_sbm(m), m)

  m$theta <- 1:10 + 0
  expect_error(
    validate_undirected_sbm(m),
    "`theta` must be equal for all nodes in a given block.",
    fixed = TRUE
  )
})

test_that("print.undirected_sbm() separates its sections", {
  set.seed(1)
  m <- sbm(n = 10, B = diag(2) * 0.5 + 0.1)
  out <- paste(capture.output(print(m)), collapse = "\n")

  expect_match(
    out,
    "Block probabilities (pi): 2 [numeric] \n\nFactor model parameterization:",
    fixed = TRUE
  )
  expect_match(out, "Poisson edges: TRUE", fixed = TRUE)
  expect_match(out, "Allow self loops: TRUE", fixed = TRUE)
})
