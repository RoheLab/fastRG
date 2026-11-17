# eigs_sym() and svds() methods agree with dense decompositions of
# expectation(), and models built from base matrices and from the various
# Matrix classes are the same model

projection <- function(U) U %*% solve(crossprod(U), t(U))

test_that("eigs_sym() matches eigen() of the expectation", {
  set.seed(60)
  n <- 30

  models <- list(
    ufm = undirected_factor_model(
      matrix(rexp(n * 3), n),
      matrix(c(0.5, 0.1, 0.2, 0.1, 0.3, 0.05, 0.2, 0.05, 0.4), 3)
    ),
    # a mixing matrix with a negative eigenvalue
    indefinite_sbm = sbm(n = n, B = matrix(c(0.1, 0.5, 0.5, 0.1), 2), block_sizes = c(12, 18)),
    dcsbm = dcsbm(theta = rexp(n) + 0.5, B = diag(3) * 0.5 + c(0.1, 0.2, 0.3)),
    mmsbm = mmsbm(theta = rexp(n) + 0.5, B = diag(3) * 0.5 + 0.1, alpha = rep(1, 3)),
    overlapping = overlapping_sbm(n = n, B = diag(3) * 0.5 + 0.1)
  )

  for (name in names(models)) {
    model <- models[[name]]
    k <- model$k
    dense <- eigen(dense_expectation(model), symmetric = TRUE)
    top <- order(abs(dense$values), decreasing = TRUE)[seq_len(k)]

    eig <- eigs_sym(model)
    expect_length(eig$values, k)

    expect_equal(
      sort(eig$values),
      sort(dense$values[top]),
      tolerance = 1e-8,
      label = name
    )
    expect_equal(
      projection(eig$vectors),
      projection(dense$vectors[, top, drop = FALSE]),
      tolerance = 1e-6,
      label = name
    )

    # fewer eigenpairs than the rank
    eig1 <- eigs_sym(model, k = 1)
    expect_equal(abs(eig1$values), max(abs(dense$values)), tolerance = 1e-8, label = name)
  }
})

test_that("svds() matches svd() of the expectation for undirected models", {
  set.seed(61)
  n <- 30
  model <- undirected_factor_model(
    matrix(rexp(n * 3), n),
    matrix(c(0.1, 0.5, 0.2, 0.5, 0.1, 0.05, 0.2, 0.05, 0.4), 3)
  )

  dense <- svd(dense_expectation(model))
  s <- svds(model)

  expect_equal(s$d, dense$d[1:3], tolerance = 1e-8)
  expect_equal(projection(s$u), projection(dense$u[, 1:3]), tolerance = 1e-6)
  expect_equal(projection(s$v), projection(dense$v[, 1:3]), tolerance = 1e-6)

  s2 <- svds(model, k = 2, nu = 0, nv = 2)
  expect_equal(s2$d, dense$d[1:2], tolerance = 1e-8)
  expect_equal(ncol(s2$v), 2)
})

test_that("svds() matches svd() of the expectation for directed models", {
  set.seed(62)

  for (d in c(30, 45)) {
    model <- directed_factor_model(
      X = matrix(rexp(30 * 3), 30),
      S = matrix(runif(6), 3, 2),
      Y = matrix(rexp(d * 2), d)
    )
    P <- dense_expectation(model)
    expect_equal(dim(P), c(30, d))

    dense <- svd(P)
    s <- svds(model)

    # the default k is the smaller rank
    expect_length(s$d, 2)
    expect_equal(s$d, dense$d[1:2], tolerance = 1e-8, label = paste("d", d))
    expect_equal(projection(s$u), projection(dense$u[, 1:2]), tolerance = 1e-6)
    expect_equal(projection(s$v), projection(dense$v[, 1:2]), tolerance = 1e-6)

    s1 <- svds(model, k = 1)
    expect_equal(s1$d, dense$d[1], tolerance = 1e-8)
  }
})

test_that("directed_dcsbm() spectra match its expectation", {
  set.seed(63)
  model <- directed_dcsbm(
    theta_out = rexp(30) + 0.5,
    theta_in = rexp(30) + 0.5,
    B = matrix(c(0.6, 0.1, 0.2, 0.1, 0.5, 0.3), 2)
  )
  dense <- svd(dense_expectation(model))
  expect_equal(svds(model)$d, dense$d[1:2], tolerance = 1e-8)
})

test_that("Matrix inputs give the same model as base matrix inputs", {
  set.seed(64)
  n <- 20
  X <- matrix(rexp(n * 2), n)
  X[sample(length(X), 15)] <- 0
  S <- matrix(c(0.4, 0.1, 0.1, 0.3), 2)

  reference <- undirected_factor_model(X, S, expected_degree = 3)

  Xs <- list(
    dgeMatrix = Matrix::Matrix(X, sparse = FALSE),
    dgCMatrix = Matrix::Matrix(X, sparse = TRUE)
  )
  Ss <- list(
    dgeMatrix = Matrix::Matrix(S, sparse = FALSE),
    dsyMatrix = Matrix::forceSymmetric(Matrix::Matrix(S, sparse = FALSE)),
    dgCMatrix = methods::as(Matrix::Matrix(S, sparse = TRUE), "generalMatrix")
  )

  for (x_name in names(Xs)) {
    for (s_name in names(Ss)) {
      label <- paste(x_name, s_name)
      m <- undirected_factor_model(Xs[[x_name]], Ss[[s_name]], expected_degree = 3)

      expect_equal(unname(as.matrix(m$S)), unname(as.matrix(reference$S)), label = label)
      expect_equal(expected_edges(m), expected_edges(reference), label = label)
      expect_equal(expected_degrees(m), expected_degrees(reference), label = label)
      expect_equal(expected_density(m), expected_density(reference), label = label)
      expect_equal(
        unname(dense_expectation(m)),
        unname(dense_expectation(reference)),
        label = label
      )

      # and the same seed gives the same graph
      set.seed(1)
      A <- sample_sparse(m)
      set.seed(1)
      A_ref <- sample_sparse(reference)
      expect_equal(A, A_ref, label = label)
    }
  }
})

test_that("diagonal Matrix mixing matrices work in every constructor", {
  D <- Matrix::Diagonal(2, c(0.4, 0.2))
  B <- as.matrix(D)

  for (poisson in c(TRUE, FALSE)) {
    set.seed(65)
    m_ddi <- sbm(n = 12, B = D, poisson_edges = poisson)
    set.seed(65)
    m_base <- sbm(n = 12, B = B, poisson_edges = poisson)
    expect_equal(unname(as.matrix(m_ddi$S)), unname(as.matrix(m_base$S)))
    expect_identical(m_ddi$z, m_base$z)

    set.seed(66)
    d_ddi <- directed_dcsbm(theta_out = rep(1, 12), theta_in = rep(1, 12), B = D, poisson_edges = poisson)
    set.seed(66)
    d_base <- directed_dcsbm(theta_out = rep(1, 12), theta_in = rep(1, 12), B = B, poisson_edges = poisson)
    expect_equal(unname(as.matrix(d_ddi$S)), unname(as.matrix(d_base$S)))
    expect_equal(expected_edges(d_ddi), expected_edges(d_base))

    set.seed(67)
    o_ddi <- overlapping_sbm(n = 12, B = D, poisson_edges = poisson)
    set.seed(67)
    o_base <- overlapping_sbm(n = 12, B = B, poisson_edges = poisson)
    expect_equal(unname(as.matrix(o_ddi$S)), unname(as.matrix(o_base$S)))
    expect_equal(unname(as.matrix(o_ddi$X)), unname(as.matrix(o_base$X)))
  }

  set.seed(68)
  m <- mmsbm(theta = rep(1, 12), B = D, alpha = c(1, 1))
  expect_equal(unname(as.matrix(m$S)), B)

  set.seed(69)
  el <- sample_edgelist(m)
  expect_s3_class(el, "tbl_df")
})

test_that("non-symmetric S is symmetrized in undirected models", {
  set.seed(70)
  X <- matrix(rexp(20), 10)
  S <- matrix(c(0.4, 0.3, 0.1, 0.2), 2)

  m <- undirected_factor_model(X, S)
  expect_equal(unname(as.matrix(m$S)), (S + t(S)) / 2)

  P <- dense_expectation(m)
  expect_equal(P, t(P))

  m2 <- dcsbm(theta = rep(1, 10), B = S, sort_nodes = FALSE)
  expect_true(isSymmetric(unname(as.matrix(m2$S))))
})
