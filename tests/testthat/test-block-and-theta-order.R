# block j of z must use row/column j of B (and pi[j]) after dcsbm() sorts
# blocks by size, and theta[i] must be the weight of node i in X

block_of <- function(X) {
  X <- as.matrix(X)
  max.col(X != 0, ties.method = "first")
}

test_that("dcsbm() keeps block labels matched to B when sorting by size", {
  B <- diag(c(0.9, 0.1, 0.1, 0.1, 0.1)) + 0.01
  sizes <- c(800, 260, 340, 280, 320)

  for (args in list(list(block_sizes = sizes), list(pi = sizes / sum(sizes)))) {
    set.seed(27)
    m <- do.call(dcsbm, c(list(theta = rep(1, 2000), B = B), args))

    S <- as.matrix(m$S)
    counts <- as.vector(table(m$z))

    # the largest block is the one the user gave the strongest
    # within-block connection, and it is still the strongest after sorting
    largest <- which.max(counts)
    expect_equal(S[largest, largest], max(diag(S)))

    # pi and block sizes line up with the relabeled z
    if (!is.null(args$block_sizes)) {
      expect_equal(counts, sort(sizes))
    }
    expect_equal(m$pi, sort(sizes) / sum(sizes))

    # z agrees with the columns of X
    expect_equal(block_of(m$X), as.integer(m$z))
  }
})

test_that("dcsbm() relabeling matches the expected edge structure", {
  set.seed(27)

  # block 1 is large with only within-block edges, block 2 is small with none
  B <- matrix(c(0.8, 0, 0, 0), 2, 2) + diag(c(0, 0.01))
  m <- dcsbm(theta = rep(1, 400), B = B, block_sizes = c(300, 100))

  A <- sample_sparse(m)
  z <- as.integer(m$z)
  big <- which.max(table(z))

  within_big <- sum(A[z == big, z == big])
  rest <- sum(A) - within_big
  expect_gt(within_big, 10 * rest)
})

test_that("dcsbm() stores theta in node order when nodes are not sorted", {
  set.seed(27)
  theta <- seq(1, 2, length.out = 60)

  for (fi in c(TRUE, FALSE)) {
    m <- dcsbm(
      theta = theta,
      B = diag(3) * 0.5 + 0.1,
      sort_nodes = FALSE,
      force_identifiability = fi
    )

    X <- as.matrix(m$X)
    expect_equal(unname(rowSums(X)), unname(m$theta))
    expect_equal(block_of(m$X), as.integer(m$z))

    if (fi) {
      # normalized within blocks, but proportional to the input in node order
      expect_equal(
        unname(m$theta),
        theta / stats::ave(theta, m$z, FUN = sum)
      )
    } else {
      expect_equal(unname(m$theta), theta)
    }
  }
})

test_that("directed_dcsbm() stores theta in node order when nodes are not sorted", {
  set.seed(27)
  theta_out <- seq(1, 2, length.out = 60)
  theta_in <- seq(3, 1, length.out = 60)

  for (fi in c(TRUE, FALSE)) {
    m <- directed_dcsbm(
      theta_out = theta_out,
      theta_in = theta_in,
      B = matrix(c(0.5, 0.1, 0.2, 0.1, 0.4, 0.3), 3, 2),
      sort_nodes = FALSE,
      force_identifiability = fi
    )

    X <- as.matrix(m$X)
    Y <- as.matrix(m$Y)

    expect_equal(unname(rowSums(X)), unname(m$theta_out))
    expect_equal(unname(rowSums(Y)), unname(m$theta_in))
    expect_equal(block_of(m$X), as.integer(m$z_out))
    expect_equal(block_of(m$Y), as.integer(m$z_in))

    if (!fi) {
      expect_equal(unname(m$theta_out), theta_out)
      expect_equal(unname(m$theta_in), theta_in)
    }
  }
})

test_that("sorted models keep X, z and theta consistent", {
  set.seed(27)
  m <- dcsbm(theta = runif(90) + 1, B = diag(3) * 0.5 + 0.1, pi = c(0.5, 0.2, 0.3))
  expect_equal(unname(rowSums(as.matrix(m$X))), unname(m$theta))
  expect_equal(block_of(m$X), as.integer(m$z))
  expect_false(is.unsorted(as.integer(m$z)))

  set.seed(27)
  d <- directed_dcsbm(
    theta_out = runif(90) + 1,
    theta_in = runif(90) + 1,
    B = matrix(runif(6), 3, 2),
    pi_out = c(0.5, 0.2, 0.3),
    pi_in = c(0.7, 0.3)
  )
  expect_equal(unname(rowSums(as.matrix(d$X))), unname(d$theta_out))
  expect_equal(unname(rowSums(as.matrix(d$Y))), unname(d$theta_in))
  expect_equal(block_of(d$X), as.integer(d$z_out))
  expect_equal(block_of(d$Y), as.integer(d$z_in))
})
