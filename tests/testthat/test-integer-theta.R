# integer theta (e.g. `rpois(n, 2) + 1L`) must give the same model as the
# equivalent double theta

test_that("dcsbm() accepts integer theta", {
  theta <- c(3L, 1L, 4L, 1L, 5L, 9L, 2L, 6L, 5L, 3L)

  for (fi in c(TRUE, FALSE)) {
    for (sort_nodes in c(TRUE, FALSE)) {
      set.seed(1)
      m_int <- dcsbm(
        theta = theta, B = diag(2), force_identifiability = fi,
        sort_nodes = sort_nodes
      )
      set.seed(1)
      m_dbl <- dcsbm(
        theta = as.numeric(theta), B = diag(2), force_identifiability = fi,
        sort_nodes = sort_nodes
      )

      expect_equal(m_int$X, m_dbl$X)
      expect_equal(m_int$theta, m_dbl$theta)
      expect_equal(m_int$z, m_dbl$z)
      expect_equal(expected_degree(m_int), expected_degree(m_dbl))
      expect_s4_class(sample_sparse(m_int), "sparseMatrix")
    }
  }
})

test_that("chung_lu() accepts integer theta", {
  theta <- c(3L, 1L, 4L, 1L, 5L, 9L, 2L, 6L, 5L, 3L)

  for (fi in c(TRUE, FALSE)) {
    for (sort_nodes in c(TRUE, FALSE)) {
      m_int <- chung_lu(
        theta = theta, force_identifiability = fi, sort_nodes = sort_nodes
      )
      m_dbl <- chung_lu(
        theta = as.numeric(theta), force_identifiability = fi,
        sort_nodes = sort_nodes
      )

      expect_equal(m_int$X, m_dbl$X)
      expect_equal(m_int$theta, m_dbl$theta)
      expect_equal(expected_degree(m_int), expected_degree(m_dbl))
      expect_s4_class(sample_sparse(m_int), "sparseMatrix")
    }
  }
})

test_that("directed_dcsbm() accepts integer theta", {
  theta_out <- c(3L, 1L, 4L, 1L, 5L, 9L, 2L, 6L, 5L, 3L)
  theta_in <- rev(theta_out)

  for (fi in c(TRUE, FALSE)) {
    for (sort_nodes in c(TRUE, FALSE)) {
      for (k in c(1, 2)) {
        set.seed(1)
        m_int <- directed_dcsbm(
          theta_out = theta_out, theta_in = theta_in, B = diag(k),
          force_identifiability = fi, sort_nodes = sort_nodes
        )
        set.seed(1)
        m_dbl <- directed_dcsbm(
          theta_out = as.numeric(theta_out), theta_in = as.numeric(theta_in),
          B = diag(k), force_identifiability = fi, sort_nodes = sort_nodes
        )

        expect_equal(m_int$X, m_dbl$X)
        expect_equal(m_int$Y, m_dbl$Y)
        expect_equal(m_int$theta_out, m_dbl$theta_out)
        expect_equal(m_int$theta_in, m_dbl$theta_in)
        expect_equal(expected_edges(m_int), expected_edges(m_dbl))
        expect_s4_class(sample_sparse(m_int), "sparseMatrix")
      }
    }
  }
})
