# internal consistency of constructed models across every combination of
# the structural flags. none of these sample graphs, so they are fast

symmetrize <- function(B) (B + t(B)) / 2

# the matrix S should be after symmetrizing, reordering and (for Bernoulli
# graphs) transforming B to the Poisson parameterization
expected_S <- function(B, o, poisson_edges) {
  S <- symmetrize(B)[o, o, drop = FALSE]
  if (poisson_edges) S else -log(1 - S)
}

test_that("dcsbm() keeps z, X, theta, pi and S consistent for every flag combination", {
  n <- 17 # not divisible by k
  theta <- seq(0.5, 3, length.out = n)

  Bs <- list(
    `1` = matrix(0.3),
    `2` = matrix(c(0.5, 0.1, 0.2, 0.3), 2), # not symmetric
    `3` = matrix(c(0.6, 0.1, 0.05, 0.1, 0.4, 0.2, 0.05, 0.2, 0.3), 3)
  )

  specs <- list(
    default = function(k) list(),
    block_sizes = function(k) list(block_sizes = c(rep(0, k - 1), n)[seq_len(k)]),
    block_sizes_mixed = function(k) list(block_sizes = balanced_block_sizes(n, k)[k:1]),
    pi = function(k) list(pi = seq_len(k) / sum(seq_len(k))),
    pi_zero = function(k) list(pi = if (k == 1) 1 else c(0, rep(1, k - 1)))
  )

  for (k_name in names(Bs)) {
    B <- Bs[[k_name]]
    k <- nrow(B)

    for (spec_name in names(specs)) {
      spec <- specs[[spec_name]](k)

      # the vector the blocks get sorted by
      sort_key <- if (!is.null(spec$block_sizes)) {
        spec$block_sizes / n
      } else if (!is.null(spec$pi)) {
        spec$pi
      } else {
        balanced_block_sizes(n, k) / n
      }

      for (sort_nodes in c(TRUE, FALSE)) {
        for (fi in c(TRUE, FALSE)) {
          for (poisson in c(TRUE, FALSE)) {
            for (loops in c(TRUE, FALSE)) {
              label <- paste(
                "k", k, spec_name, "sort", sort_nodes, "fi", fi,
                "poisson", poisson, "loops", loops
              )

              set.seed(1)
              m <- do.call(
                dcsbm,
                c(
                  list(
                    theta = theta,
                    B = B,
                    sort_nodes = sort_nodes,
                    force_identifiability = fi,
                    poisson_edges = poisson,
                    allow_self_loops = loops
                  ),
                  spec
                )
              )

              expect_s3_class(m, "undirected_dcsbm")
              expect_dcsbm_consistent(m)

              # stored fields match the arguments
              expect_identical(m$sorted, sort_nodes, label = label)
              expect_identical(m$poisson_edges, poisson, label = label)
              expect_identical(m$allow_self_loops, loops, label = label)
              expect_equal(m$n, n, label = label)
              expect_equal(m$k, k, label = label)

              o <- if (k > 1 && sort_nodes) order(sort_key) else seq_len(k)

              # B was permuted consistently with pi
              expect_equal(m$pi, sort_key[o] / sum(sort_key), label = label)
              expect_equal(
                unname(as.matrix(m$S)),
                unname(expected_S(B, o, poisson)),
                tolerance = 1e-12,
                label = label
              )

              counts <- as.vector(table(m$z))

              # exact block sizes are respected, after relabeling
              if (!is.null(spec$block_sizes)) {
                expect_equal(counts, spec$block_sizes[o], label = label)
              }

              # zero-probability blocks get no nodes
              expect_true(all(counts[m$pi == 0] == 0), label = label)

              if (fi) {
                sums <- tapply(m$theta, m$z, sum)
                expect_equal(
                  as.vector(sums[counts > 0]),
                  rep(1, sum(counts > 0)),
                  label = label
                )
              }

              if (sort_nodes) {
                expect_false(is.unsorted(as.integer(m$z)), label = label)
                # within each block, theta is non-increasing
                for (block in split(m$theta, m$z)) {
                  expect_false(is.unsorted(rev(block)), label = label)
                }
              } else if (!fi) {
                # unsorted models keep theta in node order
                expect_equal(m$theta, theta, label = label)
              }
            }
          }
        }
      }
    }
  }
})

test_that("dcsbm() with sort_nodes = FALSE keeps theta paired with its node", {
  # with force_identifiability = TRUE, theta[i] is the input theta[i]
  # divided by the total theta of node i's block
  set.seed(2)
  theta <- rexp(30) + 0.1
  m <- dcsbm(
    theta = theta,
    B = diag(3) * 0.2 + 0.05,
    sort_nodes = FALSE,
    force_identifiability = TRUE
  )

  block_totals <- tapply(theta, m$z, sum)
  expect_equal(m$theta, theta / block_totals[as.integer(m$z)], ignore_attr = TRUE)
})

test_that("sbm() and planted_partition() have theta = 1 and the right S", {
  for (sort_nodes in c(TRUE, FALSE)) {
    for (poisson in c(TRUE, FALSE)) {
      set.seed(3)
      B <- matrix(c(0.5, 0.1, 0.1, 0.3), 2)
      m <- sbm(n = 11, B = B, sort_nodes = sort_nodes, poisson_edges = poisson)
      expect_s3_class(m, "undirected_sbm")
      expect_dcsbm_consistent(m)
      expect_equal(m$theta, rep(1, 11))
      expect_equal(unname(as.matrix(expectation(m))), unname(dense_expectation(m)))

      set.seed(3)
      pp <- planted_partition(
        n = 11,
        k = 3,
        within_block = 0.4,
        between_block = 0.05,
        sort_nodes = sort_nodes,
        poisson_edges = poisson
      )
      expect_dcsbm_consistent(pp)

      S <- as.matrix(pp$S)
      transform <- if (poisson) identity else function(p) -log(1 - p)
      expect_equal(unique(diag(S)), transform(0.4))
      expect_equal(unique(S[row(S) != col(S)]), transform(0.05))
      expect_equal(sort(as.vector(table(pp$z))), c(3, 4, 4))
    }
  }
})

test_that("chung_lu() is a one-block dcsbm with theta as given", {
  theta <- c(5, 1, 3, 2, 4, 2)

  for (sort_nodes in c(TRUE, FALSE)) {
    for (fi in c(TRUE, FALSE)) {
      set.seed(4)
      m <- chung_lu(theta = theta, sort_nodes = sort_nodes, force_identifiability = fi)

      expect_s3_class(m, "undirected_dcsbm")
      expect_dcsbm_consistent(m)
      expect_equal(m$k, 1)
      expect_equal(nlevels(m$z), 1)
      expect_equal(as.matrix(m$S), matrix(1))

      expected_theta <- if (fi) theta / sum(theta) else theta
      if (sort_nodes) expected_theta <- sort(expected_theta, decreasing = TRUE)
      expect_equal(m$theta, expected_theta)

      # expectation is the outer product of theta
      expect_equal(
        unname(dense_expectation(m)),
        outer(expected_theta, expected_theta)
      )
    }
  }
})

test_that("chung_lu() generates theta from n with a message", {
  set.seed(5)
  expect_message(m <- chung_lu(n = 8), "LogNormal")
  expect_equal(m$n, 8)
  expect_length(m$theta, 8)
  expect_error(chung_lu(), "Must specify either `n` or `theta`")
  expect_error(chung_lu(n = 0), "positive integer")
})

test_that("mmsbm() keeps Z, theta, X and alpha consistent", {
  n <- 25
  theta <- seq(1, 2, length.out = n)
  B <- matrix(c(0.5, 0.1, 0.05, 0.1, 0.4, 0.2, 0.05, 0.2, 0.3), 3)
  alpha <- c(3, 1, 2)
  o <- order(alpha)

  for (sort_nodes in c(TRUE, FALSE)) {
    for (force_pure in c(TRUE, FALSE)) {
      for (poisson in c(TRUE, FALSE)) {
        label <- paste("sort", sort_nodes, "pure", force_pure, "poisson", poisson)
        set.seed(6)
        m <- mmsbm(
          theta = theta,
          B = B,
          alpha = alpha,
          sort_nodes = sort_nodes,
          force_pure = force_pure,
          poisson_edges = poisson
        )

        expect_s3_class(m, "undirected_mmsbm")
        expect_identical(m$sorted, sort_nodes)
        expect_equal(m$alpha, sort(alpha), label = label)
        expect_equal(dim(m$Z), c(n, 3), label = label)
        expect_equal(unname(rowSums(m$Z)), rep(1, n), label = label)
        expect_true(all(m$Z >= 0), label = label)

        # X = diag(theta) Z, row by row
        expect_equal(
          unname(as.matrix(m$X)),
          unname(as.matrix(m$theta * m$Z)),
          label = label
        )
        expect_setequal(m$theta, theta)
        if (!sort_nodes) expect_equal(m$theta, theta, label = label)

        expect_equal(
          unname(as.matrix(m$S)),
          unname(expected_S(B, o, poisson)),
          tolerance = 1e-12,
          label = label
        )

        if (force_pure) {
          # every block has at least one pure node
          expect_true(all(colSums(m$Z == 1) >= 1), label = label)
        }

        if (sort_nodes) {
          expect_equal(sort_by_all_columns_idx(m$X), seq_len(n), label = label)
        }
      }
    }
  }
})

test_that("mmsbm() validates its inputs", {
  expect_error(mmsbm(B = diag(2)), "Must specify either `n` or `theta`")
  expect_error(mmsbm(n = 0, B = diag(2)), "positive integer")
  expect_error(mmsbm(theta = rep(1, 10)), "Must specify either `k` or `B`")
  expect_error(mmsbm(theta = rep(1, 10), k = 0, alpha = 1), "positive integer")
  expect_error(mmsbm(theta = rep(1, 10), B = matrix(1, 2, 3)), "square")
  expect_error(
    mmsbm(theta = rep(1, 10), B = diag(2), alpha = c(1, 1, 1)),
    "must match length of `alpha`"
  )

  set.seed(7)
  expect_message(mmsbm(theta = rep(1, 10), k = 2), "Uniform")
  expect_message(mmsbm(n = 10, B = diag(2)), "LogNormal")
})

test_that("overlapping_sbm() keeps Z, X, pi and B consistent", {
  n <- 30
  B <- matrix(c(0.6, 0.1, 0.05, 0.1, 0.4, 0.2, 0.05, 0.2, 0.3), 3)
  pi <- c(0.5, 0.2, 0.3)
  o <- order(pi)

  for (sort_nodes in c(TRUE, FALSE)) {
    for (force_pure in c(TRUE, FALSE)) {
      for (poisson in c(TRUE, FALSE)) {
        for (loops in c(TRUE, FALSE)) {
          label <- paste(
            "sort", sort_nodes, "pure", force_pure,
            "poisson", poisson, "loops", loops
          )
          set.seed(8)
          m <- overlapping_sbm(
            n = n,
            B = B,
            pi = pi,
            sort_nodes = sort_nodes,
            force_pure = force_pure,
            poisson_edges = poisson,
            allow_self_loops = loops
          )

          X <- as.matrix(m$X)
          expect_identical(m$sorted, sort_nodes, label = label)
          expect_identical(m$poisson_edges, poisson, label = label)
          expect_identical(m$allow_self_loops, loops, label = label)
          expect_equal(dim(X), c(n, 3), label = label)
          expect_true(all(X %in% c(0, 1)), label = label)
          expect_true(all(rowSums(X) >= 1), label = label)
          expect_equal(unname(as.matrix(m$Z)), unname(X), label = label)
          expect_equal(m$pi, sort(pi), label = label)
          expect_equal(unname(as.matrix(m$B)), B[o, o], label = label)
          expect_equal(
            unname(as.matrix(m$S)),
            unname(expected_S(B, o, poisson)),
            tolerance = 1e-12,
            label = label
          )

          if (force_pure) {
            expect_true(all(colSums(X[rowSums(X) == 1, , drop = FALSE]) >= 1))
          }

          if (sort_nodes) {
            expect_equal(sort_by_all_columns_idx(X), seq_len(n), label = label)
          }
        }
      }
    }
  }
})

test_that("overlapping_sbm() validates its inputs", {
  expect_error(overlapping_sbm(n = 10), "Must specify either `k` or `B`")
  expect_error(overlapping_sbm(n = 10, k = 0), "positive integer")
  expect_error(overlapping_sbm(n = 10, B = matrix(0.1, 2, 3)), "square")
  expect_error(
    overlapping_sbm(n = 10, B = diag(2) * 0.5, pi = c(0.5, 0.5, 0.5)),
    "Length of `pi`"
  )
  expect_error(
    overlapping_sbm(n = 10, B = diag(2) * 0.5, pi = c(0.5, 1.5)),
    "contained in \\[0, 1\\]"
  )
  expect_error(overlapping_sbm(n = 10, B = diag(2) * 1.5), "contained in \\[0, 1\\]")
  expect_error(overlapping_sbm(n = 10, B = matrix(0.5, 2, 2)), "invertible")

  set.seed(9)
  expect_message(m <- overlapping_sbm(n = 20, k = 3), "0.8 on the diagonal")
  expect_equal(unname(diag(as.matrix(m$B))), rep(0.8, 3))
  expect_equal(unique(as.matrix(m$B)[row(m$B) != col(m$B)]), 0.05)
})

test_that("directed_dcsbm() is consistent for every flag combination", {
  n <- 13
  theta_out <- seq(0.5, 2, length.out = n)
  theta_in <- rev(seq(1, 4, length.out = n))

  Bs <- list(
    square_1 = matrix(0.4),
    square_3 = matrix(c(0.6, 0.1, 0.05, 0.2, 0.4, 0.2, 0.05, 0.3, 0.3), 3),
    rect_1_3 = matrix(c(0.2, 0.3, 0.4), 1),
    rect_3_2 = matrix(c(0.6, 0.1, 0.05, 0.2, 0.4, 0.3), 3)
  )

  for (B_name in names(Bs)) {
    B <- Bs[[B_name]]

    for (shared_blocks in c(FALSE, TRUE)) {
      if (shared_blocks && nrow(B) != ncol(B)) next

      for (shared_theta in c(FALSE, TRUE)) {
        if (shared_theta && !shared_blocks) next

        for (sort_nodes in c(TRUE, FALSE)) {
          for (fi in c(TRUE, FALSE)) {
            for (poisson in c(TRUE, FALSE)) {
              label <- paste(
                B_name, "shared_blocks", shared_blocks,
                "shared_theta", shared_theta, "sort", sort_nodes,
                "fi", fi, "poisson", poisson
              )

              theta_args <- if (shared_theta) {
                list(theta_out = theta_out)
              } else {
                list(theta_out = theta_out, theta_in = theta_in)
              }

              set.seed(10)
              m <- do.call(
                directed_dcsbm,
                c(
                  theta_args,
                  list(
                    B = B,
                    shared_blocks = shared_blocks,
                    shared_theta = shared_theta,
                    sort_nodes = sort_nodes,
                    force_identifiability = fi,
                    poisson_edges = poisson
                  )
                )
              )

              expect_s3_class(m, "directed_dcsbm")
              expect_directed_dcsbm_consistent(m)
              expect_identical(m$shared_blocks, shared_blocks, label = label)
              expect_identical(m$shared_theta, shared_theta, label = label)
              expect_identical(m$sorted, sort_nodes, label = label)
              expect_identical(m$poisson_edges, poisson, label = label)

              # default block sizes are balanced on both sides
              expect_equal(
                sort(as.vector(table(m$z_out))),
                sort(balanced_block_sizes(n, nrow(B))),
                label = label
              )
              expect_equal(
                sort(as.vector(table(m$z_in))),
                sort(balanced_block_sizes(n, ncol(B))),
                label = label
              )

              if (shared_blocks) {
                expect_identical(m$z_out, m$z_in, label = label)
              }

              if (shared_theta) {
                expect_identical(m$theta_out, m$theta_in, label = label)
              }

              if (fi) {
                expect_equal(
                  as.vector(tapply(m$theta_out, m$z_out, sum)),
                  rep(1, m$k1),
                  label = label
                )
                expect_equal(
                  as.vector(tapply(m$theta_in, m$z_in, sum)),
                  rep(1, m$k2),
                  label = label
                )
              }

              if (sort_nodes) {
                expect_false(is.unsorted(as.integer(m$z_out)), label = label)
                if (!shared_blocks) {
                  expect_false(is.unsorted(as.integer(m$z_in)), label = label)
                }
              } else if (!fi) {
                expect_equal(m$theta_out, theta_out, label = label)
                expected_in <- if (shared_theta) theta_out else theta_in
                expect_equal(m$theta_in, expected_in, label = label)
              }

              transform <- if (poisson) identity else function(p) -log(1 - p)
              expect_setequal(as.vector(as.matrix(m$S)), transform(as.vector(B)))
            }
          }
        }
      }
    }
  }
})

test_that("directed_dcsbm() k_out / k_in generate a random B", {
  set.seed(11)
  expect_message(
    m <- directed_dcsbm(theta_out = rep(1, 10), theta_in = rep(1, 10), k_out = 2, k_in = 3),
    "Uniform"
  )
  expect_equal(dim(m$S), c(2, 3))

  # with shared blocks, either k_out or k_in is enough
  set.seed(11)
  m <- suppressMessages(
    directed_dcsbm(theta_out = rep(1, 10), theta_in = rep(1, 10), k_out = 2, shared_blocks = TRUE)
  )
  expect_equal(dim(m$S), c(2, 2))
})

test_that("directed_dcsbm() validates theta, k and pi", {
  B <- diag(2) * 0.5
  th <- rep(1, 10)

  expect_error(
    directed_dcsbm(n = 0, B = B, shared_blocks = TRUE, shared_theta = TRUE),
    "positive integer"
  )
  expect_error(
    directed_dcsbm(n = 10, theta_out = rep(1, 10), theta_in = rep(1, 10), B = B),
    "not all three at once"
  )
  expect_error(directed_dcsbm(n = 0, B = B), "positive integer")
  expect_error(
    directed_dcsbm(theta_out = rep(1, 10), theta_in = rep(1, 9), B = B),
    "Length of `theta_out` must match length of `theta_in`"
  )
  expect_error(
    directed_dcsbm(n = 10, theta_out = rep(1, 10), B = B),
    "but not both"
  )
  expect_error(directed_dcsbm(theta_out = th, theta_in = th), "Must specify either `k_out` and `k_in` together")
  expect_error(directed_dcsbm(theta_out = th, theta_in = th, k_out = 2), "Must specify either `k_out` and `k_in` together")
  expect_error(directed_dcsbm(theta_out = th, theta_in = th, k_out = 0, k_in = 2), "`k_out` must be a positive integer")
  expect_error(directed_dcsbm(theta_out = th, theta_in = th, k_out = 2, k_in = 0), "`k_in` must be a positive integer")
  expect_error(
    directed_dcsbm(theta_out = th, theta_in = th, B = B, pi_out = c(1, 1, 1)),
    "Length of `pi_out` must match"
  )
  expect_error(
    directed_dcsbm(theta_out = th, theta_in = th, B = B, pi_in = c(1, 1, 1)),
    "Length of `pi_in` must match"
  )
  expect_error(
    directed_dcsbm(theta_out = th, theta_in = th, B = B, pi_out = c(-1, 2)),
    "`pi_out` must be >= 0"
  )
  expect_error(
    directed_dcsbm(theta_out = th, theta_in = th, B = B, pi_in = c(-1, 2)),
    "`pi_in` must be >= 0"
  )

  set.seed(12)
  expect_message(m <- directed_dcsbm(n = 10, B = B), "LogNormal")
  expect_length(m$theta_out, 10)
})

test_that("random theta and B come with messages in dcsbm()", {
  set.seed(13)
  expect_message(m <- dcsbm(n = 10, B = diag(2)), "LogNormal")
  expect_length(m$theta, 10)
  expect_message(m <- dcsbm(theta = rep(1, 10), k = 3), "Uniform")
  expect_equal(m$k, 3)
})

test_that("dcsbm() validates its inputs", {
  expect_error(dcsbm(B = diag(2)), "Must specify either `n` or `theta`")
  expect_error(dcsbm(n = 0, B = diag(2)), "positive integer")
  expect_error(dcsbm(n = 10, theta = rep(1, 10), B = diag(2)), "only one of `n` and `theta`")
  expect_error(dcsbm(theta = rep(1, 10)), "Must specify either `k` or `B`")
  expect_error(dcsbm(theta = rep(1, 10), k = 0), "positive integer")
  expect_error(dcsbm(theta = rep(1, 10), B = matrix(1, 2, 3)), "square")
  expect_error(
    dcsbm(theta = rep(1, 10), B = diag(2), block_sizes = c(5, 6)),
    "Sum of `block_sizes`"
  )
  expect_error(
    dcsbm(theta = rep(1, 10), B = diag(2), pi = c(1, 1, 1)),
    "Length of `pi`"
  )
  expect_error(
    dcsbm(theta = rep(1, 10), B = diag(2), pi = c(-1, 2)),
    "must be >= 0"
  )
  expect_error(dcsbm(theta = c(-1, rep(1, 9)), B = diag(2)), "non-negative")
  expect_error(dcsbm(theta = rep(1, 10), B = -diag(2)), "non-negative")
})
