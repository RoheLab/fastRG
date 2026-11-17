empty_rectangular <- function() {
  directed_factor_model(
    X = matrix(1, 20, 2),
    S = diag(2),
    Y = matrix(1, 15, 2),
    expected_density = 1e-9
  )
}

test_that("sample_sparse() handles empty rectangular directed graphs", {
  set.seed(27)
  A <- sample_sparse(empty_rectangular())

  expect_s4_class(A, "dgCMatrix")
  expect_equal(dim(A), c(20, 15))
  expect_equal(Matrix::nnzero(A), 0)
})

test_that("sample_igraph() and sample_tidygraph() handle empty rectangular graphs", {
  set.seed(27)
  g <- sample_igraph(empty_rectangular())
  expect_equal(igraph::vcount(g), 35)
  expect_equal(igraph::ecount(g), 0)

  set.seed(27)
  tg <- sample_tidygraph(empty_rectangular())
  expect_equal(igraph::vcount(tg), 35)
  expect_equal(igraph::ecount(tg), 0)
})

test_that("sample_sparse() returns the same class for empty and non-empty graphs", {
  X <- matrix(1, 20, 2)

  set.seed(27)
  empty_dir <- sample_sparse(
    directed_factor_model(X, diag(2), X, expected_density = 1e-9)
  )
  full_dir <- sample_sparse(
    directed_factor_model(X, diag(2), X, expected_density = 0.5)
  )
  expect_equal(Matrix::nnzero(empty_dir), 0)
  expect_equal(class(empty_dir), class(full_dir))

  empty_undir <- sample_sparse(
    undirected_factor_model(X, diag(2), expected_density = 1e-9)
  )
  full_undir <- sample_sparse(
    undirected_factor_model(X, diag(2), expected_density = 0.5)
  )
  expect_equal(Matrix::nnzero(empty_undir), 0)
  expect_equal(dim(empty_undir), c(20, 20))
  expect_equal(class(empty_undir), class(full_undir))
})

test_that("sample_edgelist() returns the same type for empty and non-empty graphs", {
  X <- matrix(1, 20, 2)
  Y <- matrix(1, 15, 2)

  models <- list(
    undirected = function(d) undirected_factor_model(X, diag(2), expected_density = d),
    directed = function(d) directed_factor_model(X, diag(2), X, expected_density = d),
    rectangular = function(d) directed_factor_model(X, diag(2), Y, expected_density = d)
  )

  for (model in models) {
    set.seed(27)
    empty <- sample_edgelist(model(1e-9))
    full <- sample_edgelist(model(0.5))

    expect_equal(nrow(empty), 0)
    expect_gt(nrow(full), 0)
    expect_s3_class(empty, "tbl_df")
    expect_identical(class(empty), class(full))
    expect_identical(vapply(empty, typeof, ""), vapply(full, typeof, ""))
    expect_identical(vapply(empty, typeof, ""), c(from = "integer", to = "integer"))
  }
})
