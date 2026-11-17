# tiny graphs, empty graphs, empty blocks and other boundary cases

test_that("tiny undirected and directed models sample with every sampler", {
  set.seed(100)
  models <- list(
    er_1 = erdos_renyi(n = 1, p = 0.5),
    er_2 = erdos_renyi(n = 2, p = 0.5),
    sbm_2 = sbm(n = 2, B = diag(2) * 0.5),
    pp_3 = planted_partition(n = 3, k = 3, within_block = 0.5, between_block = 0.1),
    directed_er_1 = directed_erdos_renyi(n = 1, p = 0.5),
    directed_2 = directed_dcsbm(theta_out = c(1, 2), theta_in = c(2, 1), B = diag(2)),
    rectangular_1_by_3 = directed_factor_model(matrix(1), matrix(0.5), matrix(1, 3, 1))
  )

  for (name in names(models)) {
    model <- models[[name]]
    directed <- inherits(model, "directed_factor_model")
    d <- if (directed) model$d else model$n

    for (seed in 1:5) {
      s <- sample_all(model, seed)
      expect_equal(dim(s$A), c(model$n, d), label = name)
      expect_s3_class(s$el, "tbl_df")
      expect_equal(igraph::ecount(s$ig), nrow(s$el), label = name)
    }
  }
})

test_that("models with no edges give empty graphs of the right shape", {
  models <- list(
    undirected = undirected_factor_model(matrix(1, 6, 1), matrix(0)),
    directed_square = directed_erdos_renyi(n = 6, p = 0),
    directed_rectangular = directed_factor_model(matrix(1, 6, 1), matrix(0), matrix(1, 4, 1)),
    zero_block_sbm = sbm(n = 6, B = matrix(c(0, 0, 0, 0), 2), poisson_edges = FALSE)
  )

  for (name in names(models)) {
    model <- models[[name]]
    expect_equal(expected_edges(model), 0, label = name)

    for (multigraph in c(TRUE, FALSE)) {
      s <- sample_all(model, 1, multigraph = multigraph)
      expect_equal(nrow(s$el), 0, label = name)
      expect_equal(sum(s$A), 0, label = name)
      expect_equal(igraph::ecount(s$ig), 0, label = name)
      expect_equal(igraph::ecount(s$tg), 0, label = name)
      if (!multigraph) {
        expect_named(s$el, c("from", "to", "weight"))
      }
    }
  }
})

test_that("empty blocks and zero-probability blocks get no nodes and no edges", {
  set.seed(101)
  B <- matrix(c(0.5, 0.2, 0.1, 0.2, 0.4, 0.3, 0.1, 0.3, 0.6), 3)

  m <- sbm(n = 20, B = B, block_sizes = c(8, 0, 12))
  expect_equal(sort(as.vector(table(m$z))), c(0, 8, 12))
  empty <- which(table(m$z) == 0)
  expect_equal(sum(as.matrix(m$X)[, empty]), 0)

  m <- dcsbm(theta = rexp(20) + 0.5, B = B, pi = c(0.5, 0, 0.5))
  empty <- which(m$pi == 0)
  expect_length(empty, 1)
  expect_equal(sum(m$z == levels(m$z)[empty]), 0)

  # sampling never touches the empty block
  P <- dense_expectation(m)
  for (seed in 1:5) {
    set.seed(seed)
    A <- sample_sparse(m)
    expect_true(all(as.matrix(A)[P == 0] == 0))
  }

  # directed models with empty blocks on either side
  set.seed(102)
  dm <- directed_dcsbm(
    theta_out = rep(1, 15),
    theta_in = rep(1, 15),
    B = matrix(runif(6), 2),
    block_sizes_out = c(0, 15),
    pi_in = c(0.5, 0, 0.5)
  )
  expect_true(any(table(dm$z_out) == 0))
  expect_true(any(table(dm$z_in) == 0))
  expect_equal(dim(sample_sparse(dm)), c(15, 15))
})

test_that("balanced default block sizes differ by at most one", {
  for (n in c(5, 7, 10, 23)) {
    for (k in 1:5) {
      sizes <- balanced_block_sizes(n, k)
      expect_equal(sum(sizes), n)
      expect_length(sizes, k)
      expect_lte(max(sizes) - min(sizes), 1)

      set.seed(103)
      m <- sbm(n = n, B = diag(k) * 0.5 + 0.1)
      expect_equal(sort(as.vector(table(m$z))), sort(sizes))
    }
  }
})

test_that("rank one models of every kind work", {
  set.seed(104)
  models <- list(
    sbm = sbm(n = 10, B = matrix(0.3)),
    dcsbm = dcsbm(theta = rexp(10), B = matrix(0.3)),
    dcsbm_unsorted = dcsbm(theta = rexp(10), B = matrix(0.3), sort_nodes = FALSE),
    planted_partition = planted_partition(n = 10, k = 1, within_block = 0.3, between_block = 0.1),
    overlapping = overlapping_sbm(n = 10, B = matrix(0.3)),
    directed = directed_dcsbm(theta_out = rexp(10), theta_in = rexp(10), B = matrix(0.3)),
    directed_pp = directed_planted_partition(n = 10, k = 1, within_block = 0.3, between_block = 0.1)
  )

  for (name in names(models)) {
    model <- models[[name]]
    P <- dense_expectation(model)
    expect_equal(qr(P)$rank, 1, label = name)
    expect_no_error(sample_sparse(model))
  }
})

test_that("planted_partition() validates its inputs", {
  expect_error(
    planted_partition(n = 2, k = 3, within_block = 0.5, between_block = 0.1),
    "less than or equal to `n`"
  )
  expect_error(
    planted_partition(n = 10, k = 2, within_block = 0.5, between_block = c(0.1, 0.2)),
    "`between_block` must be a scalar"
  )
  expect_error(
    planted_partition(n = 10, k = 2, within_block = 0.5, between_block = 1.5),
    "`between_block` must be between zero and one"
  )
  expect_error(
    planted_partition(n = 10, k = 2, within_block = -0.5, between_block = 0.1),
    "non-negative"
  )
})

test_that("planted_partition() a and b parameterization divides by n", {
  set.seed(105)
  m <- planted_partition(n = 50, k = 2, a = 10, b = 2)
  expect_equal(m$within_block, 10 / 50)
  expect_equal(m$between_block, 2 / 50)
  S <- as.matrix(m$S)
  expect_equal(unique(diag(S)), 0.2)
  expect_equal(unique(S[row(S) != col(S)]), 0.04)
})

test_that("sample_edgelist() works directly on matrices", {
  set.seed(106)
  X <- matrix(rexp(10), 5)
  S <- matrix(c(0.5, 0.1, 0.2, 0.3), 2)
  Y <- matrix(rexp(14), 7)

  el <- sample_edgelist(X, S, Y, TRUE, TRUE, TRUE)
  expect_s3_class(el, "tbl_df")
  expect_true(all(el$from %in% 1:5))
  expect_true(all(el$to %in% 1:7))

  set.seed(106)
  el_matrix <- sample_edgelist(Matrix::Matrix(X), Matrix::Matrix(S), Matrix::Matrix(Y), TRUE, TRUE, TRUE)
  set.seed(106)
  el_base <- sample_edgelist(X, S, Y, TRUE, TRUE, TRUE)
  expect_identical(el_matrix, el_base)

  # undirected sampling puts every edge in the upper triangle
  el <- sample_edgelist(X, S, X, FALSE, TRUE, TRUE)
  expect_true(all(el$from <= el$to))

  expect_error(sample_edgelist(X, S, Y, "yes", TRUE, TRUE))
  expect_error(sample_edgelist(X, S, Y, TRUE, TRUE, TRUE, multigraph = NA), "multigraph")
})
