test_that("plot_sparse_matrix() has the same orientation as plot_dense_matrix()", {
  # a single entry in row 1, column 3 of a 2 x 4 matrix
  A <- Matrix::sparseMatrix(i = 1, j = 3, x = 1, dims = c(2, 4))

  sparse <- ggplot2::layer_data(plot_sparse_matrix(A))
  dense <- ggplot2::layer_data(plot_dense_matrix(as.matrix(A)))
  dense <- dense[dense$fill == max(dense$fill), ]

  # columns along the x axis, rows down the (reversed) y axis
  expect_equal(sparse$x, 3)
  expect_equal(sparse$y, -1)
  expect_equal(sparse$x, dense$x)
  expect_equal(sparse$y, dense$y)
})
