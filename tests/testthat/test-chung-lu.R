test_that("chung_lu() returns an undirected_chung_lu object", {
  set.seed(6)
  m <- chung_lu(theta = rexp(15) + 0.5)

  expect_s3_class(
    m,
    c("undirected_chung_lu", "undirected_dcsbm", "undirected_factor_model"),
    exact = TRUE
  )

  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(out, "Undirected Chung-Lu Graph", fixed = TRUE)
  expect_match(out, "Degree heterogeneity (theta): 15 [numeric]", fixed = TRUE)
  expect_no_match(out, "Degree-Corrected", fixed = TRUE)
})
