# theta = 0 is allowed and gives isolated nodes. sorting must keep each
# node's block, theta and row of X together

check_alignment <- function(X, z, theta, sorted) {
  X <- as.matrix(X)
  z <- as.integer(z)

  expect_equal(unname(rowSums(X)), unname(theta))

  # every row has at most one non-zero entry, in the column of its block
  nonzero <- theta > 0
  expect_equal(
    max.col(X != 0, ties.method = "first")[nonzero],
    z[nonzero]
  )
  expect_true(all(rowSums(X != 0) <= 1))

  if (sorted) {
    expect_false(is.unsorted(z))
    for (b in unique(z)) {
      expect_false(is.unsorted(-theta[z == b]))
    }
  }
}

test_that("dcsbm() keeps z, X and theta aligned with zero theta", {
  theta <- c(0, 1, 2, 3, 0, 1, 2, 3, 4, 5)

  for (sort_nodes in c(TRUE, FALSE)) {
    for (fi in c(TRUE, FALSE)) {
      for (k in 1:2) {
        set.seed(3)
        m <- dcsbm(
          theta = theta, B = diag(k), sort_nodes = sort_nodes,
          force_identifiability = fi
        )
        check_alignment(m$X, m$z, m$theta, sort_nodes)
        expect_false(anyNA(m$theta))

        # the zero-theta nodes are the ones with zero expected degree
        expect_equal(sum(m$theta == 0), 2)
        expect_equal(
          which(m$theta == 0),
          unname(which(Matrix::rowSums(m$X %*% m$S %*% t(m$X)) == 0))
        )

        A <- sample_sparse(m)
        expect_equal(sum(A[m$theta == 0, ]), 0)
      }
    }
  }
})

test_that("directed_dcsbm() keeps z, X/Y and theta aligned with zero theta", {
  theta_out <- c(0, 1, 2, 3, 0, 1, 2, 3, 4, 5)
  theta_in <- c(5, 0, 4, 3, 2, 0, 1, 3, 2, 1)

  for (sort_nodes in c(TRUE, FALSE)) {
    for (fi in c(TRUE, FALSE)) {
      for (shared_blocks in c(TRUE, FALSE)) {
        set.seed(3)
        m <- directed_dcsbm(
          theta_out = theta_out, theta_in = theta_in, B = diag(2),
          sort_nodes = sort_nodes, force_identifiability = fi,
          shared_blocks = shared_blocks
        )
        check_alignment(m$X, m$z_out, m$theta_out, sort_nodes && !shared_blocks)
        check_alignment(m$Y, m$z_in, m$theta_in, sort_nodes && !shared_blocks)

        A <- sample_sparse(m)
        expect_equal(sum(A[m$theta_out == 0, ]), 0)
        expect_equal(sum(A[, m$theta_in == 0]), 0)
      }
    }
  }
})

test_that("negative theta is still an error", {
  expect_error(
    dcsbm(theta = c(-1, 1, 1, 1), B = diag(2)),
    "non-negative"
  )
})
