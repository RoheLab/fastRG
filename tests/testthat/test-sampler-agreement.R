# sample_edgelist(), sample_sparse(), sample_igraph() and sample_tidygraph()
# describe the same graph for a common seed, and every sampled graph obeys
# the structural constraints of its model

agreement_models <- function() {
  set.seed(40)
  n <- 30
  list(
    ufm = undirected_factor_model(
      matrix(rexp(n * 2), n),
      matrix(c(0.2, 0.05, 0.1, 0.3), 2),
      expected_degree = 6
    ),
    erdos_renyi = erdos_renyi(n = n, p = 0.2),
    sbm_k1 = sbm(n = n, B = matrix(0.3)),
    dcsbm = dcsbm(theta = rexp(n) + 0.5, B = diag(3) * 0.4 + 0.05, expected_degree = 5),
    dcsbm_unsorted = dcsbm(
      theta = rexp(n) + 0.5,
      B = diag(2) * 0.4 + 0.05,
      sort_nodes = FALSE,
      expected_degree = 5
    ),
    mmsbm = mmsbm(theta = rexp(n) + 0.5, B = diag(2) * 0.4 + 0.05, alpha = c(1, 1)),
    overlapping = overlapping_sbm(n = n, B = diag(2) * 0.4 + 0.05),
    directed_square = directed_factor_model(
      matrix(rexp(n * 2), n),
      matrix(runif(6), 2, 3),
      matrix(rexp(n * 3), n),
      expected_out_degree = 5
    ),
    directed_rect_wide = directed_factor_model(
      matrix(rexp(n * 2), n),
      matrix(runif(4), 2),
      matrix(rexp(2 * (n + 7)), n + 7),
      expected_out_degree = 5
    ),
    directed_rect_tall = directed_factor_model(
      matrix(rexp(n), n),
      matrix(1),
      matrix(rexp(11), 11),
      expected_in_degree = 4
    ),
    directed_dcsbm = directed_dcsbm(
      theta_out = rexp(n),
      theta_in = rexp(n),
      B = matrix(runif(6), 2),
      expected_out_degree = 5
    ),
    directed_pp = directed_planted_partition(
      n = n,
      k = 3,
      within_block = 0.4,
      between_block = 0.05
    )
  )
}

test_that("all samplers agree and respect allow_self_loops / poisson_edges", {
  models <- agreement_models()

  for (name in names(models)) {
    for (loops in c(TRUE, FALSE)) {
      for (poisson in c(TRUE, FALSE)) {
        model <- models[[name]]
        model$allow_self_loops <- loops
        model$poisson_edges <- poisson
        label <- paste(name, "loops", loops, "poisson", poisson)

        directed <- inherits(model, "directed_factor_model")
        rectangular <- directed && model$n != model$d
        n <- model$n
        d <- if (directed) model$d else n

        for (seed in 1:3) {
          s <- sample_all(model, seed)
          el <- s$el
          A <- s$A
          ig <- s$ig
          tg <- s$tg

          # output types
          expect_s3_class(el, "tbl_df")
          expect_type(el$from, "integer")
          expect_type(el$to, "integer")
          expect_equal(dim(A), c(n, d), label = label)
          if (directed) {
            expect_s4_class(A, "dgCMatrix")
          } else {
            expect_s4_class(A, "dsCMatrix")
          }
          expect_s3_class(tg, "tbl_graph")

          # node indices are in range
          expect_true(all(el$from >= 1 & el$from <= n), label = label)
          expect_true(all(el$to >= 1 & el$to <= d), label = label)

          # the edgelist and the adjacency matrix agree
          diag_A <- if (rectangular) 0 else sum(Matrix::diag(A))
          if (directed) {
            expect_equal(nrow(el), sum(A), label = label)
          } else {
            expect_true(all(el$from <= el$to), label = label)
            expect_equal(nrow(el), (sum(A) + diag_A) / 2, label = label)
          }

          # igraph and tidygraph have every node and every edge
          expect_equal(igraph::ecount(ig), nrow(el), label = label)
          expect_equal(igraph::ecount(tg), nrow(el), label = label)
          expected_vertices <- if (rectangular) n + d else n
          expect_equal(igraph::vcount(ig), expected_vertices, label = label)
          expect_equal(igraph::vcount(tg), expected_vertices, label = label)

          # degrees agree with A
          if (rectangular) {
            deg <- igraph::degree(ig)
            expect_equal(unname(deg[seq_len(n)]), unname(rowSums(A)), label = label)
            expect_equal(unname(deg[n + seq_len(d)]), unname(colSums(A)), label = label)
          } else if (directed) {
            expect_equal(igraph::is_directed(ig), TRUE)
            expect_equal(unname(igraph::degree(ig, mode = "out")), unname(rowSums(A)), label = label)
            expect_equal(unname(igraph::degree(ig, mode = "in")), unname(colSums(A)), label = label)
          } else {
            expect_false(igraph::is_directed(ig))
            # igraph counts each self-loop twice towards the degree
            expect_equal(
              unname(igraph::degree(ig)),
              unname(rowSums(A) + Matrix::diag(A)),
              label = label
            )
          }

          # structural constraints
          if (!loops) {
            expect_false(any(el$from == el$to), label = label)
            if (!rectangular) {
              expect_equal(sum(Matrix::diag(A)), 0, label = label)
            }
          }

          if (!poisson) {
            expect_false(anyDuplicated(el[c("from", "to")]) > 0, label = label)
            expect_lte(max(c(0, A@x)), 1)
            expect_false(any(igraph::which_multiple(ig)), label = label)
          }

          # collapsed weights add up to the multigraph edge count
          set.seed(seed)
          weighted <- sample_edgelist(model, multigraph = FALSE)
          expect_equal(sum(weighted$weight), nrow(el), label = label)
          expect_false(anyDuplicated(weighted[c("from", "to")]) > 0, label = label)
          if (!poisson) {
            expect_true(all(weighted$weight == 1), label = label)
          }
        }
      }
    }
  }
})

test_that("sampling is reproducible for a fixed seed", {
  models <- agreement_models()

  for (name in names(models)) {
    model <- models[[name]]

    set.seed(99)
    el1 <- sample_edgelist(model)
    set.seed(99)
    el2 <- sample_edgelist(model)
    expect_identical(el1, el2, label = name)

    set.seed(99)
    A1 <- sample_sparse(model)
    set.seed(99)
    A2 <- sample_sparse(model)
    expect_identical(A1, A2, label = name)
  }
})

test_that("constructors are reproducible for a fixed seed", {
  builders <- list(
    dcsbm = function() dcsbm(n = 20, k = 3),
    sbm = function() sbm(n = 20, k = 2),
    chung_lu = function() chung_lu(n = 20),
    mmsbm = function() mmsbm(n = 20, k = 3),
    overlapping = function() overlapping_sbm(n = 20, k = 3),
    directed = function() directed_dcsbm(n = 20, k_out = 2, k_in = 3),
    directed_pp = function() {
      directed_planted_partition(n = 20, k = 2, a = 4, b = 1, pi = c(0.3, 0.7))
    }
  )

  for (name in names(builders)) {
    set.seed(100)
    m1 <- suppressMessages(builders[[name]]())
    set.seed(100)
    m2 <- suppressMessages(builders[[name]]())
    expect_identical(m1, m2, label = name)
  }
})

test_that("sorting nodes only relabels the nodes of an sbm", {
  # with theta = 1, the sorted and unsorted models from the same seed have
  # the same blocks, so their expectations agree up to a permutation
  B <- matrix(c(0.5, 0.1, 0.05, 0.1, 0.4, 0.2, 0.05, 0.2, 0.3), 3)

  for (args in list(list(pi = c(0.5, 0.2, 0.3)), list(block_sizes = c(9, 4, 7)))) {
    set.seed(41)
    sorted <- do.call(sbm, c(list(n = 20, B = B, sort_nodes = TRUE), args))
    set.seed(41)
    unsorted <- do.call(sbm, c(list(n = 20, B = B, sort_nodes = FALSE), args))

    P_sorted <- dense_expectation(sorted)
    P_unsorted <- dense_expectation(unsorted)

    # match blocks by their within-block rate, which is distinct per block
    block_rate_sorted <- diag(as.matrix(sorted$S))[as.integer(sorted$z)]
    block_rate_unsorted <- diag(as.matrix(unsorted$S))[as.integer(unsorted$z)]
    expect_equal(sort(block_rate_sorted), sort(block_rate_unsorted))

    perm <- order(-block_rate_unsorted)
    perm_sorted <- order(-block_rate_sorted)
    expect_equal(
      unname(P_sorted[perm_sorted, perm_sorted]),
      unname(P_unsorted[perm, perm])
    )
  }
})
