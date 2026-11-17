# all samplers consume the RNG identically (sample_sparse(), sample_igraph()
# and sample_tidygraph() all go through sample_edgelist()), so with a common
# seed the collapsed weights must match sample_sparse() entry by entry

# nonzero entries of a sparse adjacency matrix as a sorted (from, to, weight)
# data frame. for undirected graphs keep only the upper triangle, matching
# the from <= to convention of sample_edgelist()
sparse_triplets <- function(A, undirected) {
  A <- as(as(A, "generalMatrix"), "TsparseMatrix")
  df <- data.frame(from = A@i + 1L, to = A@j + 1L, weight = A@x)
  df <- df[df$weight != 0, ]
  if (undirected) {
    df <- df[df$from <= df$to, ]
  }
  sort_triplets(df)
}

sort_triplets <- function(df) {
  df <- as.data.frame(df)
  df$from <- as.integer(df$from)
  df$to <- as.integer(df$to)
  df$weight <- as.integer(df$weight)
  df <- df[order(df$from, df$to), c("from", "to", "weight")]
  rownames(df) <- NULL
  df
}

# edges of a one-mode igraph as a sorted (from, to, weight) data frame
igraph_triplets <- function(g, undirected) {
  df <- igraph::as_data_frame(g, what = "edges")
  from <- as.integer(df$from)
  to <- as.integer(df$to)
  if (undirected) {
    df$from <- pmin(from, to)
    df$to <- pmax(from, to)
  } else {
    df$from <- from
    df$to <- to
  }
  sort_triplets(df)
}

# edges of a bipartite igraph built by graph_from_biadjacency_matrix(), where
# vertices 1..n are rows of A and vertices n + 1, ..., n + d are columns
bipartite_triplets <- function(g, n) {
  df <- igraph::as_data_frame(g, what = "edges")
  from <- pmin(df$from, df$to)
  to <- pmax(df$from, df$to)
  df$from <- from
  df$to <- to - n
  sort_triplets(df)
}

make_ufm <- function(poisson_edges = TRUE, allow_self_loops = TRUE) {
  set.seed(30)
  n <- 60
  k <- 3
  X <- matrix(rpois(n = n * k, 1), nrow = n)
  S <- matrix(runif(n = k * k, 0, .1), nrow = k)
  # dense enough that there are many repeated edges and self-loops
  undirected_factor_model(
    X,
    S,
    expected_degree = 20,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )
}

make_fm <- function(d = 60, poisson_edges = TRUE, allow_self_loops = TRUE) {
  set.seed(31)
  n <- 60
  k1 <- 3
  k2 <- 2
  X <- matrix(rpois(n = n * k1, 1), nrow = n)
  S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
  Y <- matrix(rexp(n = k2 * d, 1), nrow = d)
  directed_factor_model(
    X,
    S,
    Y,
    expected_in_degree = 20,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops
  )
}

test_that("multigraph = TRUE is the default and is unchanged", {
  ufm <- make_ufm()
  fm_square <- make_fm(d = 60)
  fm_rect <- make_fm(d = 40)

  for (model in list(ufm, fm_square, fm_rect)) {
    set.seed(1)
    el_default <- sample_edgelist(model)
    set.seed(1)
    el_true <- sample_edgelist(model, multigraph = TRUE)
    expect_identical(el_default, el_true)
    expect_named(el_default, c("from", "to"))

    # the multigraph really does have repeated edges
    expect_true(anyDuplicated(el_default) > 0)

    set.seed(1)
    ig_default <- sample_igraph(model)
    set.seed(1)
    ig_true <- sample_igraph(model, multigraph = TRUE)
    expect_true(igraph::identical_graphs(ig_default, ig_true))
    expect_false("weight" %in% igraph::edge_attr_names(ig_default))
    expect_true(any(igraph::which_multiple(ig_default)))

    # one edge per sampled edge
    set.seed(1)
    A <- sample_sparse(model)
    total <- if (inherits(model, "undirected_factor_model")) {
      sum(Matrix::triu(A))
    } else {
      sum(A)
    }
    expect_equal(igraph::ecount(ig_default), total)
    expect_equal(nrow(el_default), total)

    set.seed(1)
    tg_default <- sample_tidygraph(model)
    expect_s3_class(tg_default, "tbl_graph")
    expect_false("weight" %in% igraph::edge_attr_names(tg_default))
    expect_equal(igraph::ecount(tg_default), total)
  }
})

test_that("undirected: collapsed weights match sample_sparse()", {
  ufm <- make_ufm()

  set.seed(2)
  A <- sample_sparse(ufm)
  expected <- sparse_triplets(A, undirected = TRUE)

  # make sure the test actually exercises collapsing and self-loops
  expect_true(any(expected$weight > 1))
  expect_true(any(expected$from == expected$to))

  set.seed(2)
  el <- sample_edgelist(ufm, multigraph = FALSE)
  expect_named(el, c("from", "to", "weight"))
  expect_type(el$weight, "integer")
  expect_true(all(el$from <= el$to))
  expect_equal(anyDuplicated(el[c("from", "to")]), 0)
  expect_equal(sort_triplets(el), expected)

  set.seed(2)
  ig <- sample_igraph(ufm, multigraph = FALSE)
  expect_false(igraph::is_directed(ig))
  expect_equal(igraph::vcount(ig), ufm$n)
  expect_type(igraph::E(ig)$weight, "integer")
  expect_false(any(igraph::which_multiple(ig)))
  expect_equal(igraph_triplets(ig, undirected = TRUE), expected)

  set.seed(2)
  tg <- sample_tidygraph(ufm, multigraph = FALSE)
  expect_s3_class(tg, "tbl_graph")
  expect_equal(igraph_triplets(tg, undirected = TRUE), expected)

  # total weight equals the number of rows of the multigraph edgelist
  set.seed(2)
  el_multi <- sample_edgelist(ufm)
  expect_equal(sum(el$weight), nrow(el_multi))
})

test_that("undirected: (i, j) and (j, i) are collapsed together", {
  ufm <- make_ufm()

  set.seed(3)
  el <- sample_edgelist(ufm, multigraph = FALSE)
  unordered <- paste(pmin(el$from, el$to), pmax(el$from, el$to))
  expect_equal(anyDuplicated(unordered), 0)
})

test_that("undirected without self-loops: weights match sample_sparse()", {
  ufm <- make_ufm(allow_self_loops = FALSE)

  set.seed(4)
  expected <- sparse_triplets(sample_sparse(ufm), undirected = TRUE)

  set.seed(4)
  el <- sample_edgelist(ufm, multigraph = FALSE)
  expect_false(any(el$from == el$to))
  expect_equal(sort_triplets(el), expected)

  set.seed(4)
  ig <- sample_igraph(ufm, multigraph = FALSE)
  expect_equal(igraph_triplets(ig, undirected = TRUE), expected)
})

test_that("directed square: collapsed weights match sample_sparse()", {
  fm <- make_fm(d = 60)

  set.seed(5)
  A <- sample_sparse(fm)
  expected <- sparse_triplets(A, undirected = FALSE)
  expect_true(any(expected$weight > 1))

  set.seed(5)
  el <- sample_edgelist(fm, multigraph = FALSE)
  expect_named(el, c("from", "to", "weight"))
  expect_type(el$weight, "integer")
  expect_equal(sort_triplets(el), expected)

  # directed edges are not symmetrized: (i, j) and (j, i) stay distinct
  expect_true(any(el$from > el$to))

  set.seed(5)
  ig <- sample_igraph(fm, multigraph = FALSE)
  expect_true(igraph::is_directed(ig))
  expect_type(igraph::E(ig)$weight, "integer")
  expect_equal(igraph_triplets(ig, undirected = FALSE), expected)

  set.seed(5)
  tg <- sample_tidygraph(fm, multigraph = FALSE)
  expect_s3_class(tg, "tbl_graph")
  expect_equal(igraph_triplets(tg, undirected = FALSE), expected)
})

test_that("directed rectangular: collapsed weights match sample_sparse()", {
  fm <- make_fm(d = 40)

  set.seed(6)
  A <- sample_sparse(fm)
  expected <- sparse_triplets(A, undirected = FALSE)
  expect_true(any(expected$weight > 1))

  set.seed(6)
  el <- sample_edgelist(fm, multigraph = FALSE)
  expect_equal(sort_triplets(el), expected)

  set.seed(6)
  ig <- sample_igraph(fm, multigraph = FALSE)
  expect_true(igraph::is_bipartite(ig))
  expect_equal(igraph::vcount(ig), fm$n + fm$d)
  expect_type(igraph::E(ig)$weight, "integer")
  expect_false(any(igraph::which_multiple(ig)))
  expect_equal(bipartite_triplets(ig, fm$n), expected)

  set.seed(6)
  tg <- sample_tidygraph(fm, multigraph = FALSE)
  expect_s3_class(tg, "tbl_graph")
  expect_equal(bipartite_triplets(tg, fm$n), expected)
})

test_that("poisson_edges = FALSE gives a weight column of 1s", {
  models <- list(
    make_ufm(poisson_edges = FALSE),
    make_fm(d = 60, poisson_edges = FALSE),
    make_fm(d = 40, poisson_edges = FALSE)
  )

  for (model in models) {
    # multigraph = TRUE: no duplicates and no weight, as before
    set.seed(7)
    el_multi <- sample_edgelist(model)
    expect_named(el_multi, c("from", "to"))
    expect_equal(anyDuplicated(el_multi), 0)

    set.seed(7)
    el <- sample_edgelist(model, multigraph = FALSE)
    expect_named(el, c("from", "to", "weight"))
    expect_true(all(el$weight == 1L))
    expect_equal(nrow(el), nrow(el_multi))

    set.seed(7)
    ig <- sample_igraph(model, multigraph = FALSE)
    expect_true(all(igraph::E(ig)$weight == 1L))
    expect_equal(igraph::ecount(ig), nrow(el_multi))

    set.seed(7)
    tg <- sample_tidygraph(model, multigraph = FALSE)
    expect_true(all(igraph::E(tg)$weight == 1L))
  }
})

test_that("multigraph = FALSE with no sampled edges", {
  X <- matrix(1, nrow = 30, ncol = 1)
  ufm <- undirected_factor_model(X, matrix(1), expected_degree = 1e-10)
  fm <- directed_factor_model(X, matrix(1), X, expected_in_degree = 1e-10)

  set.seed(8)
  el <- sample_edgelist(ufm, multigraph = FALSE)
  expect_equal(nrow(el), 0)
  expect_named(el, c("from", "to", "weight"))

  set.seed(8)
  el <- sample_edgelist(fm, multigraph = FALSE)
  expect_equal(nrow(el), 0)
  expect_named(el, c("from", "to", "weight"))

  set.seed(8)
  ig <- sample_igraph(fm, multigraph = FALSE)
  expect_equal(igraph::ecount(ig), 0)
  expect_equal(igraph::vcount(ig), 30)
})

test_that("multigraph is passed through and does not trip check_dots", {
  ufm <- make_ufm()

  # multigraph must be named, and unknown named arguments still error
  expect_error(sample_edgelist(ufm, foo = 1))
  expect_error(sample_igraph(ufm, foo = 1))
  expect_error(sample_tidygraph(ufm, foo = 1))

  # sample_sparse() does not take a multigraph argument
  expect_error(sample_sparse(ufm, multigraph = FALSE))
})

test_that("bad multigraph input gives a clear error", {
  ufm <- make_ufm()
  fm <- make_fm(d = 40)

  bad <- list(NA, NULL, "no", 1, c(TRUE, FALSE), logical(0))

  for (value in bad) {
    expect_error(
      sample_edgelist(ufm, multigraph = value),
      "`multigraph` must be a single `TRUE` or `FALSE`"
    )
    expect_error(
      sample_igraph(fm, multigraph = value),
      "`multigraph` must be a single `TRUE` or `FALSE`"
    )
    expect_error(
      sample_tidygraph(ufm, multigraph = value),
      "`multigraph` must be a single `TRUE` or `FALSE`"
    )
  }

  expect_snapshot(sample_igraph(fm, multigraph = NA), error = TRUE)
})
