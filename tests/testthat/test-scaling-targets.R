# expected_degree / expected_density / expected_in_degree /
# expected_out_degree rescale S to hit their target exactly, for every
# constructor, and invalid targets are rejected

undirected_builders <- function(n, ...) {
  theta <- seq(0.5, 2, length.out = n)
  list(
    undirected_factor_model = function(...) {
      undirected_factor_model(cbind(theta, rev(theta)), matrix(c(0.3, 0.1, 0.2, 0.4), 2), ...)
    },
    erdos_renyi = function(...) erdos_renyi(n = n, ...),
    sbm = function(...) sbm(n = n, B = diag(3) * 0.3 + 0.1, ...),
    dcsbm = function(...) dcsbm(theta = theta, B = diag(3) * 0.3 + 0.1, ...),
    planted_partition = function(...) {
      planted_partition(n = n, k = 3, within_block = 0.3, between_block = 0.1, ...)
    },
    chung_lu = function(...) chung_lu(theta = theta, ...),
    mmsbm = function(...) mmsbm(theta = theta, B = diag(3) * 0.3 + 0.1, alpha = rep(1, 3), ...),
    overlapping_sbm = function(...) overlapping_sbm(n = n, B = diag(3) * 0.3 + 0.1, ...)
  )
}

test_that("expected_degree and expected_density hit their targets for undirected models", {
  n <- 40
  builders <- undirected_builders(n)

  for (name in names(builders)) {
    for (loops in c(TRUE, FALSE)) {
      set.seed(80)
      m <- builders[[name]](expected_degree = 3.5, allow_self_loops = loops)
      expect_equal(expected_degree(m), 3.5, label = name)
      expect_equal(mean(expected_degrees(m)), 3.5, label = name)
      expect_identical(m$allow_self_loops, loops)

      set.seed(80)
      m <- builders[[name]](expected_density = 0.07, allow_self_loops = loops)
      expect_equal(expected_density(m), 0.07, label = name)
    }
  }
})

test_that("scaling Bernoulli blockmodels hits the target on the probability scale", {
  # for blockmodels with theta = 1, the Bernoulli edge probabilities are
  # 1 - exp(-S), and these are what expected_degree / expected_density
  # rescale
  n <- 30

  for (scaling in list(list(expected_degree = 4), list(expected_density = 0.1))) {
    builders <- list(
      sbm = function(...) sbm(n = n, B = diag(3) * 0.3 + 0.1, block_sizes = c(5, 10, 15), ...),
      planted_partition = function(...) {
        planted_partition(n = n, k = 3, within_block = 0.3, between_block = 0.1, ...)
      },
      erdos_renyi = function(...) erdos_renyi(n = n, ...)
    )

    for (name in names(builders)) {
      set.seed(81)
      m <- do.call(builders[[name]], c(scaling, poisson_edges = FALSE))
      X <- as.matrix(m$X)
      P_bernoulli <- X %*% (1 - exp(-as.matrix(m$S))) %*% t(X)

      if (!is.null(scaling$expected_degree)) {
        expect_equal(mean(rowSums(P_bernoulli)), 4, label = name)
      } else {
        off <- P_bernoulli[row(P_bernoulli) != col(P_bernoulli)]
        expect_equal(mean(off), 0.1, label = name)
      }

      # the Poisson-scale summaries overestimate slightly, as documented
      if (!is.null(scaling$expected_degree)) {
        expect_gt(expected_degree(m), 4)
      } else {
        expect_gt(expected_density(m), 0.1)
      }
    }
  }
})

test_that("directed scaling hits its target for every constructor, square or not", {
  set.seed(82)
  n <- 30
  d <- 45

  builders <- list(
    factor_model_square = function(...) {
      directed_factor_model(matrix(rexp(n * 2), n), matrix(runif(6), 2, 3), matrix(rexp(n * 3), n), ...)
    },
    factor_model_rectangular = function(...) {
      directed_factor_model(matrix(rexp(n * 2), n), matrix(runif(6), 2, 3), matrix(rexp(d * 3), d), ...)
    },
    erdos_renyi = function(...) directed_erdos_renyi(n = n, ...),
    directed_dcsbm = function(...) {
      directed_dcsbm(theta_out = rexp(n), theta_in = rexp(n), B = matrix(runif(6), 2), ...)
    },
    directed_dcsbm_shared = function(...) {
      directed_dcsbm(
        theta_out = rexp(n),
        B = matrix(runif(4), 2),
        shared_blocks = TRUE,
        shared_theta = TRUE,
        ...
      )
    },
    directed_pp = function(...) {
      directed_planted_partition(n = n, k = 3, within_block = 0.4, between_block = 0.1, ...)
    }
  )

  for (name in names(builders)) {
    m <- builders[[name]](expected_in_degree = 2.5)
    expect_equal(expected_in_degree(m), 2.5, label = name)
    P <- dense_expectation(m)
    expect_equal(mean(colSums(P)), 2.5, label = name)

    m <- builders[[name]](expected_out_degree = 3.5)
    expect_equal(expected_out_degree(m), 3.5, label = name)
    P <- dense_expectation(m)
    expect_equal(mean(rowSums(P)), 3.5, label = name)

    m <- builders[[name]](expected_density = 0.06)
    expect_equal(expected_density(m), 0.06, label = name)
    P <- dense_expectation(m)
    expect_equal(mean(P), 0.06, label = name)

    # in and out degree are tied together by the number of rows and columns
    expect_equal(
      expected_in_degree(m) * m$d,
      expected_out_degree(m) * m$n,
      label = name
    )
  }
})

test_that("undirected_factor_model() rejects invalid scaling targets", {
  X <- matrix(1, 5, 1)
  S <- matrix(0.5)

  expect_error(
    undirected_factor_model(X, S, expected_degree = 1, expected_density = 0.1),
    "Cannot specify both"
  )
  expect_error(undirected_factor_model(X, S, expected_degree = 0), "strictly greater than zero")
  expect_error(undirected_factor_model(X, S, expected_degree = -1), "strictly greater than zero")
  expect_error(undirected_factor_model(X, S, expected_density = 0), "strictly between zero and one")
  expect_error(undirected_factor_model(X, S, expected_density = 1), "strictly between zero and one")
  expect_error(
    undirected_factor_model(X, S, expected_degree = 10, poisson_edges = FALSE),
    "must not exceed 1 for Bernoulli graphs"
  )
  expect_error(undirected_factor_model(-X, S), "non-negative")
  expect_error(undirected_factor_model(X, -S), "non-negative")
  expect_error(undirected_factor_model(matrix(1, 5, 2), S), "number of columns of `X`")
  expect_error(
    undirected_factor_model(X, S, allow_self_loops = "yes"),
    "`allow_self_loops` must be a logical"
  )
})

test_that("directed_factor_model() rejects invalid scaling targets", {
  X <- matrix(1, 5, 1)
  S <- matrix(0.5)
  Y <- matrix(1, 4, 1)

  expect_error(
    directed_factor_model(X, S, Y, expected_in_degree = 1, expected_out_degree = 1),
    "Must specify only one"
  )
  expect_error(
    directed_factor_model(X, S, Y, expected_in_degree = 1, expected_density = 0.1),
    "Must specify only one"
  )
  expect_error(directed_factor_model(X, S, Y, expected_in_degree = 0), "strictly greater than zero")
  expect_error(directed_factor_model(X, S, Y, expected_out_degree = 0), "strictly greater than zero")
  expect_error(directed_factor_model(X, S, Y, expected_density = 0), "strictly between zero and one")
  expect_error(directed_factor_model(X, S, Y, expected_density = 1), "strictly between zero and one")
  expect_error(
    directed_factor_model(X, S, Y, expected_out_degree = 5, poisson_edges = FALSE),
    "must not exceed 1 for Bernoulli graphs"
  )
  expect_error(directed_factor_model(-X, S, Y), "non-negative")
  expect_error(directed_factor_model(X, S, -Y), "non-negative")
  expect_error(directed_factor_model(matrix(1, 5, 2), S, Y), "`k1` must equal")
  expect_error(directed_factor_model(X, S, matrix(1, 4, 2)), "`k2` must equal")
  expect_error(
    directed_factor_model(X, S, Y, allow_self_loops = "yes"),
    "`allow_self_loops` must be a logical"
  )
})

test_that("erdos_renyi() models have the requested edge probability", {
  for (poisson in c(TRUE, FALSE)) {
    m <- erdos_renyi(n = 25, p = 0.3, poisson_edges = poisson)
    expect_s3_class(m, "undirected_erdos_renyi")
    expect_equal(m$p, 0.3)
    S <- as.matrix(m$S)
    expect_equal(if (poisson) S else 1 - exp(-S), matrix(0.3))
    expect_equal(dim(m$X), c(25, 1))

    dm <- directed_erdos_renyi(n = 25, p = 0.3, poisson_edges = poisson)
    expect_s3_class(dm, "directed_erdos_renyi")
    expect_equal(dm$p, 0.3)
    S <- as.matrix(dm$S)
    expect_equal(if (poisson) S else 1 - exp(-S), matrix(0.3))
    expect_equal(c(dm$n, dm$d), c(25, 25))
  }

  expect_error(erdos_renyi(n = 10, p = 0), "strictly between zero and one")
  expect_error(erdos_renyi(n = 10, p = 1), "strictly between zero and one")
  expect_error(directed_erdos_renyi(n = 10, p = -0.1), "non-negative")
})
