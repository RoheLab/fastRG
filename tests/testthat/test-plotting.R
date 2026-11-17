test_that("plot_dense_matrix() draws one tile per entry with the entry's value", {
  A <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2) # 2 x 3, not symmetric

  p <- plot_dense_matrix(A)
  expect_s3_class(p, "ggplot")

  long <- p$data
  expect_equal(nrow(long), length(A))
  expect_equal(long$value, A[cbind(long$row, long$col)])

  # columns of A run along the x axis, rows down the (reversed) y axis
  built <- ggplot2::layer_data(p)
  expect_equal(nrow(built), length(A))
  expect_equal(sort(unique(built$x)), 1:3)
  expect_equal(sort(unique(-built$y)), 1:2)
})

test_that("plot_dense_matrix() ignores existing dimnames", {
  A <- matrix(1:4, 2, dimnames = list(c("a", "b"), c("c", "d")))
  p <- plot_dense_matrix(A)
  expect_equal(p$data$value, c(1, 2, 3, 4))
  expect_type(p$data$row, "double")
  expect_type(p$data$col, "double")
})

test_that("plot_expectation() plots expectation() for every kind of model", {
  set.seed(90)
  models <- list(
    undirected = sbm(n = 6, B = matrix(c(0.5, 0.1, 0.1, 0.3), 2)),
    directed_square = directed_erdos_renyi(n = 5, p = 0.2),
    directed_rectangular = directed_factor_model(
      matrix(rexp(8), 4),
      matrix(runif(4), 2),
      matrix(rexp(14), 7)
    )
  )

  for (name in names(models)) {
    model <- models[[name]]
    P <- dense_expectation(model)
    p <- plot_expectation(model)
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), length(P), label = name)
    expect_equal(p$data$value, P[cbind(p$data$row, p$data$col)], label = name)
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("plot_sparse_matrix() draws one tile per nonzero entry", {
  set.seed(91)
  undirected <- sample_sparse(sbm(n = 12, B = matrix(c(0.6, 0.1, 0.1, 0.4), 2)))
  directed <- sample_sparse(
    directed_factor_model(matrix(rexp(4), 4), matrix(2), matrix(rexp(7), 7))
  )

  for (A in list(undirected, directed)) {
    p <- plot_sparse_matrix(A)
    expect_s3_class(p, "ggplot")

    general <- as(as(A, "generalMatrix"), "TsparseMatrix")
    nonzero <- general@x != 0
    built <- ggplot2::layer_data(p)
    expect_equal(nrow(built), sum(nonzero))
    expect_setequal(as.numeric(as.character(p$data$x)), unique(general@x[nonzero]))
  }

  expect_error(plot_sparse_matrix(matrix(1, 2, 2)))
})
