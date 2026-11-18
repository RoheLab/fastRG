# models with a single node, or a single block, must keep matrix-valued
# X, Y and Z

test_that("dcsbm() and sbm() work with a single node", {
  m <- dcsbm(theta = 1, B = matrix(1))
  expect_s4_class(m$X, "Matrix")
  expect_equal(dim(m$X), c(1, 1))
  expect_equal(expected_edges(m), 1)
  expect_equal(dim(sample_sparse(m)), c(1, 1))

  m <- dcsbm(theta = 2, B = matrix(1), sort_nodes = FALSE)
  expect_equal(unname(as.matrix(m$X)), matrix(2))

  set.seed(1)
  m <- sbm(n = 1, B = diag(2))
  expect_s4_class(m$X, "Matrix")
  expect_equal(dim(m$X), c(1, 2))
  expect_equal(dim(sample_sparse(m)), c(1, 1))

  m <- erdos_renyi(n = 1, p = 0.5)
  expect_equal(dim(m$X), c(1, 1))
})

test_that("directed models work with a single node", {
  m <- directed_dcsbm(
    theta_out = 1, theta_in = 1, B = matrix(1),
    force_identifiability = FALSE
  )
  expect_equal(dim(m$X), c(1, 1))
  expect_equal(dim(m$Y), c(1, 1))
  expect_equal(dim(sample_sparse(m)), c(1, 1))

  m <- directed_erdos_renyi(n = 1, p = 0.5)
  expect_equal(dim(m$X), c(1, 1))
})

test_that("overlapping_sbm() keeps Z a matrix with a single block", {
  set.seed(1)
  m <- overlapping_sbm(n = 10, B = matrix(0.5), pi = 1)
  expect_s4_class(m$Z, "Matrix")
  expect_equal(dim(m$Z), c(10, 1))
  expect_equal(dim(m$X), c(10, 1))
})
