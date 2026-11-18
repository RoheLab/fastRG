# eigs_sym() and svds() methods for factor models use RSpectra, which is
# unreliable on small, exactly low-rank operators
# (https://github.com/yixuan/RSpectra/issues/20). these tests cover the
# smaller default ncv, the small-model warning, and argument forwarding

dense_top_singular <- function(model, k) {
  svd(as.matrix(expectation(model)))$d[seq_len(k)]
}

dense_top_abs_eigen <- function(model, k) {
  values <- eigen(as.matrix(expectation(model)), symmetric = TRUE)$values
  sort(abs(values), decreasing = TRUE)[seq_len(k)]
}

test_that("small models give correct spectra with the default ncv", {
  # with RSpectra's default ncv, these were wrong or errored in most seeds:
  # directed svds was wrong in 30 of 40 seeds, undirected svds errored in 33
  for (seed in 1:15) {
    set.seed(seed)
    n <- 20

    directed <- directed_factor_model(
      matrix(runif(2 * n), n),
      matrix(runif(4), 2),
      matrix(runif(2 * n), n)
    )
    undirected <- undirected_factor_model(
      matrix(runif(3 * n), n),
      crossprod(matrix(runif(9), 3))
    )

    expect_warning(s_dir <- svds(directed, k = 2), class = "fastRG_warning_small_model_spectra")
    expect_equal(s_dir$d, dense_top_singular(directed, 2), tolerance = 1e-6, label = seed)

    expect_warning(s_und <- svds(undirected, k = 3), class = "fastRG_warning_small_model_spectra")
    expect_equal(s_und$d, dense_top_singular(undirected, 3), tolerance = 1e-6, label = seed)

    expect_warning(e_und <- eigs_sym(undirected, k = 3), class = "fastRG_warning_small_model_spectra")
    expect_equal(
      sort(abs(e_und$values), decreasing = TRUE),
      dense_top_abs_eigen(undirected, 3),
      tolerance = 1e-6,
      label = seed
    )
  }
})

test_that("small models warn and link to the RSpectra issue", {
  set.seed(27)
  small <- sbm(n = 25, B = diag(2) * 0.5 + 0.1)

  expect_warning(eigs_sym(small), "github.com/yixuan/RSpectra/issues/20", fixed = TRUE)
  expect_warning(svds(small), class = "fastRG_warning_small_model_spectra")

  # the warning depends on the smaller dimension of directed models
  rectangular <- directed_factor_model(
    matrix(runif(400), 200),
    diag(2),
    matrix(runif(40), 20)
  )
  expect_warning(svds(rectangular), class = "fastRG_warning_small_model_spectra")
})

test_that("larger models don't warn", {
  set.seed(27)
  undirected <- sbm(n = 26, B = diag(2) * 0.5 + 0.1)
  directed <- directed_factor_model(
    matrix(runif(200), 100),
    diag(2),
    matrix(runif(60), 30)
  )

  expect_no_warning(eigs_sym(undirected))
  expect_no_warning(svds(undirected))
  expect_no_warning(svds(directed))
})

test_that("eigs_sym() forwards `which`", {
  set.seed(27)
  # B has a negative eigenvalue, so the smallest eigenvalue of E[A] is negative
  model <- sbm(n = 40, B = matrix(c(0.1, 0.5, 0.5, 0.1), 2), block_sizes = c(20, 20))
  values <- eigen(as.matrix(expectation(model)), symmetric = TRUE)$values

  expect_equal(eigs_sym(model, k = 1, which = "SA")$values, min(values), tolerance = 1e-8)
  expect_equal(eigs_sym(model, k = 1, which = "LA")$values, max(values), tolerance = 1e-8)
  expect_equal(
    eigs_sym(model, k = 2, which = "LM")$values,
    values[order(abs(values), decreasing = TRUE)][1:2],
    tolerance = 1e-8
  )
})

test_that("opts are forwarded and a user-supplied ncv is respected", {
  set.seed(27)
  model <- sbm(n = 40, B = diag(3) * 0.5 + 0.1)

  # retvec = FALSE drops the eigenvectors, which only happens if opts arrive
  expect_null(eigs_sym(model, opts = list(retvec = FALSE))$vectors)

  expect_equal(default_ncv(list(), k = 3, dim = 40)$ncv, 10)
  expect_equal(default_ncv(list(), k = 6, dim = 40)$ncv, 13)
  expect_equal(default_ncv(list(), k = 2, dim = 8)$ncv, 7)
  expect_equal(default_ncv(list(ncv = 25), k = 3, dim = 40)$ncv, 25)

  # when k < ncv <= dim is impossible, leave ncv to RSpectra
  expect_null(default_ncv(list(), k = 2, dim = 3)$ncv)
})

test_that("sigma gives an informative error for factor models", {
  model <- sbm(n = 40, B = diag(2) * 0.5 + 0.1)
  expect_error(eigs_sym(model, sigma = 1), "shift-and-invert")
})
