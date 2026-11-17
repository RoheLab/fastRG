# every print method runs and reports values that agree with the model

expect_prints_undirected_summary <- function(m) {
  out <- paste(capture.output(print(m)), collapse = "\n")

  expect_match(out, paste0("Nodes \\(n\\): ", m$n), label = class(m)[1])
  expect_match(out, paste0("Expected edges: ", round(expected_edges(m))), fixed = TRUE)
  expect_match(
    out,
    paste0("Expected degree: ", round(expected_degree(m), 1)),
    fixed = TRUE
  )
  expect_match(
    out,
    paste0("Expected density: ", round(expected_density(m), 5)),
    fixed = TRUE
  )
  expect_match(out, paste0("Poisson edges: ", m$poisson_edges), fixed = TRUE)
  expect_match(out, paste0("Allow self loops: ", m$allow_self_loops), fixed = TRUE)
  expect_match(out, paste0("X: ", m$n, " x ", ncol(m$X)), fixed = TRUE)
  expect_match(out, paste0("S: ", m$k, " x ", m$k), fixed = TRUE)

  invisible(out)
}

expect_prints_directed_summary <- function(m) {
  out <- paste(capture.output(print(m)), collapse = "\n")

  expect_match(out, paste0("Expected edges: ", round(expected_edges(m))), fixed = TRUE)
  expect_match(
    out,
    paste0("Expected in degree: ", round(expected_in_degree(m), 1)),
    fixed = TRUE
  )
  expect_match(
    out,
    paste0("Expected out degree: ", round(expected_out_degree(m), 1)),
    fixed = TRUE
  )
  expect_match(
    out,
    paste0("Expected density: ", round(expected_density(m), 5)),
    fixed = TRUE
  )
  expect_match(out, paste0("X: ", m$n, " x ", m$k1), fixed = TRUE)
  expect_match(out, paste0("S: ", m$k1, " x ", m$k2), fixed = TRUE)
  expect_match(out, paste0("Y: ", m$d, " x ", m$k2), fixed = TRUE)
  expect_match(out, paste0("Poisson edges: ", m$poisson_edges), fixed = TRUE)
  expect_match(out, paste0("Allow self loops: ", m$allow_self_loops), fixed = TRUE)

  invisible(out)
}

test_that("print.undirected_factor_model() reports the model", {
  set.seed(1)
  m <- undirected_factor_model(
    matrix(rexp(40), 20),
    diag(2) * 0.3 + 0.1,
    poisson_edges = FALSE,
    allow_self_loops = FALSE
  )
  out <- expect_prints_undirected_summary(m)
  expect_match(out, "Undirected Factor Model")
  expect_match(out, "Rank (k): 2", fixed = TRUE)

  # erdos_renyi() uses the generic undirected print method
  out <- expect_prints_undirected_summary(erdos_renyi(n = 30, p = 0.2))
  expect_match(out, "Undirected Factor Model")
})

test_that("print.directed_factor_model() reports the model", {
  set.seed(2)
  m <- directed_factor_model(
    X = matrix(rexp(30), 15),
    S = matrix(runif(6), 2, 3),
    Y = matrix(rexp(24), 8),
    expected_density = 0.1
  )
  out <- expect_prints_directed_summary(m)
  expect_match(out, "Directed Factor Model")
  expect_match(out, "(n): 15", fixed = TRUE)
  expect_match(out, "(d): 8", fixed = TRUE)
  expect_match(out, "(k1): 2", fixed = TRUE)
  expect_match(out, "(k2): 3", fixed = TRUE)

  out <- expect_prints_directed_summary(directed_erdos_renyi(n = 12, p = 0.3))
  expect_match(out, "Directed Factor Model")
})

test_that("print.undirected_dcsbm() reports sorting and parameters", {
  for (sort_nodes in c(TRUE, FALSE)) {
    set.seed(3)
    m <- dcsbm(theta = rexp(25) + 0.5, B = diag(3) * 0.4 + 0.05, sort_nodes = sort_nodes)
    out <- expect_prints_undirected_summary(m)
    expect_match(out, "Undirected Degree-Corrected Stochastic Blockmodel")
    expect_match(out, "Blocks (k): 3", fixed = TRUE)
    expect_match(out, "Block memberships (z): 25 [factor]", fixed = TRUE)
    expect_match(out, "Degree heterogeneity (theta): 25 [numeric]", fixed = TRUE)
    expect_match(out, "Block probabilities (pi): 3 [numeric]", fixed = TRUE)

    sorted <- if (sort_nodes) "(arranged by block)" else "(not arranged by block)"
    expect_match(out, paste("Nodes (n): 25", sorted), fixed = TRUE)
  }
})

test_that("print.undirected_sbm() reports sorting and parameters", {
  for (sort_nodes in c(TRUE, FALSE)) {
    set.seed(4)
    m <- sbm(n = 20, B = diag(2) * 0.4 + 0.05, sort_nodes = sort_nodes)
    out <- expect_prints_undirected_summary(m)
    expect_match(out, "Undirected Stochastic Blockmodel")
    expect_match(out, "Blocks (k): 2", fixed = TRUE)
    expect_match(out, "Block memberships (z): 20 [factor]", fixed = TRUE)

    sorted <- if (sort_nodes) "(arranged by block)" else "(not arranged by block)"
    expect_match(out, paste("Nodes (n): 20", sorted), fixed = TRUE)
  }
})

test_that("print.undirected_mmsbm() reports sorting and parameters", {
  for (sort_nodes in c(TRUE, FALSE)) {
    set.seed(5)
    m <- mmsbm(
      theta = rexp(20) + 0.5,
      B = diag(3) * 0.4 + 0.05,
      alpha = c(1, 2, 3),
      sort_nodes = sort_nodes
    )
    out <- expect_prints_undirected_summary(m)
    expect_match(out, "Mixed Membership Stochastic Blockmodel")
    expect_match(out, "Blocks (k): 3", fixed = TRUE)
    expect_match(out, "(Z): 20 x 3", fixed = TRUE)
    expect_match(out, "Degree heterogeneity (theta): 20 [numeric]", fixed = TRUE)
    expect_match(out, "Block propensities (alpha): 3 [numeric]", fixed = TRUE)

    sorted <- if (sort_nodes) "(arranged by block)" else "(not arranged by block)"
    expect_match(out, paste("Nodes (n): 20", sorted), fixed = TRUE)
  }
})

test_that("print.undirected_chung_lu() runs and reports the model", {
  # chung_lu() objects currently dispatch to print.undirected_dcsbm(), so
  # call the Chung-Lu method directly to check that it works
  set.seed(6)
  m <- chung_lu(theta = rexp(15) + 0.5)
  expect_prints_undirected_summary(m)

  out <- paste(capture.output(print.undirected_chung_lu(m)), collapse = "\n")
  expect_match(out, "Undirected Chung-Lu Graph")
  expect_match(out, "Nodes (n): 15 (arranged by block)", fixed = TRUE)
  expect_match(out, "Degree heterogeneity (theta): 15 [numeric]", fixed = TRUE)
  expect_match(out, paste0("Expected edges: ", round(expected_edges(m))), fixed = TRUE)
  expect_match(
    out,
    paste0("Expected density: ", round(expected_density(m), 5)),
    fixed = TRUE
  )
})

test_that("print methods of the remaining models report the right values", {
  set.seed(7)
  expect_prints_undirected_summary(
    planted_partition(n = 20, k = 2, within_block = 0.3, between_block = 0.05)
  )
  expect_prints_undirected_summary(
    overlapping_sbm(n = 20, B = diag(2) * 0.4 + 0.05, poisson_edges = FALSE)
  )
  expect_prints_directed_summary(
    directed_dcsbm(theta_out = rexp(20), theta_in = rexp(20), B = matrix(runif(6), 2))
  )
  expect_prints_directed_summary(
    directed_planted_partition(n = 20, k = 2, within_block = 0.3, between_block = 0.05)
  )
})

test_that("dim_and_class() describes matrices and vectors", {
  expect_equal(dim_and_class(matrix(0, 2, 3)), "2 x 3 [matrix]")
  expect_equal(dim_and_class(Matrix::Diagonal(2)), "2 x 2 [ddiMatrix]")
  expect_equal(dim_and_class(1:4), "4 [integer]")
  expect_equal(dim_and_class(factor(c("a", "b"))), "2 [factor]")
})
