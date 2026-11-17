# Verification of expected_edges(), expected_degree() and expected_density()
# for undirected factor models.
#
# Write P = E[A | X, S] = X S X'. The argument has three steps, each tested
# separately below:
#
# 1. Algebra. The closed forms used by the package equal the definitions
#    written as explicit sums over the entries of a dense P:
#
#      expected_edges   = sum_{i <= j} P[i, j]           = (sum(P) + tr(P)) / 2
#      expected_degree  = mean_i sum_j P[i, j]           = sum(P) / n
#      expected_density = mean_{i != j} P[i, j]          = (sum(P) - tr(P)) / (n (n - 1))
#
# 2. Counting. In every sampled undirected graph, each edge between
#    distinct nodes appears twice in the symmetric adjacency matrix A and
#    each self-loop appears once, so
#
#      number of edges            = (sum(A) + tr(A)) / 2
#      edges between distinct nodes = (sum(A) - tr(A)) / 2
#
#    This is a deterministic identity, checked on every sample.
#
# 3. Distribution. sample_edgelist() draws m ~ Poisson(sum(X (S / 2) X'))
#    and places each edge independently with probability proportional to
#    P, so by Poisson splitting each unordered pair {i, j}, i < j, gets an
#    independent Poisson(P[i, j] / 2 + P[j, i] / 2) = Poisson(P[i, j])
#    number of edges. Self-loops get Poisson(P[i, i] / 2) from the sampler
#    plus Poisson(P[i, i] / 2) from the top-up added in #46, which is
#    Poisson(P[i, i]). Hence, for Poisson graphs, *exactly*:
#
#      number of edges              ~ Poisson(expected_edges)
#      edges between distinct nodes ~ Poisson(expected_density * choose(n, 2))
#
#    so totals over many replicates can be tested with exact z-scores,
#    rather than ad hoc tolerances.

dense_reference <- function(model) {
  P <- as.matrix(model$X %*% model$S %*% Matrix::t(model$X))
  n <- nrow(P)
  off_diagonal <- row(P) != col(P)

  list(
    edges = sum(P[upper.tri(P, diag = TRUE)]),
    degree = mean(rowSums(P)),
    degrees = rowSums(P),
    density = mean(P[off_diagonal]),
    off_diagonal_edges = sum(P[upper.tri(P)])
  )
}

# z-score of an observed total of `reps` independent Poisson(mu) draws
poisson_z <- function(total, mu, reps) {
  (total - reps * mu) / sqrt(reps * mu)
}

reference_models <- function() {
  set.seed(27)
  n <- 60
  k <- 3

  X_dense <- matrix(rexp(n * k), nrow = n)
  S_asym <- matrix(runif(k * k), nrow = k)

  X_sparse <- Matrix::rsparsematrix(n, k, density = 0.3, rand.x = function(m) runif(m))
  X_sparse[cbind(seq_len(n), rep_len(seq_len(k), n))] <- 1

  list(
    dense_asymmetric_S = undirected_factor_model(X_dense, S_asym),
    sparse_X = undirected_factor_model(X_sparse, diag(k) * 0.3 + 0.05),
    rank_one = undirected_factor_model(matrix(rexp(n), ncol = 1), matrix(0.2)),
    scaled_degree = undirected_factor_model(X_dense, S_asym, expected_degree = 4),
    scaled_density = undirected_factor_model(X_dense, S_asym, expected_density = 0.1),
    erdos_renyi = erdos_renyi(n = n, p = 0.1),
    sbm = sbm(n = n, B = matrix(c(0.3, 0.05, 0.05, 0.2), 2), block_sizes = c(25, 35)),
    dcsbm = dcsbm(theta = rexp(n) + 0.5, B = diag(3) * 0.4 + 0.05, pi = c(0.2, 0.3, 0.5)),
    planted_partition = planted_partition(n = n, k = 3, within_block = 0.3, between_block = 0.05),
    chung_lu = chung_lu(theta = rexp(n) + 1),
    mmsbm = suppressMessages(mmsbm(theta = rexp(n) + 0.5, B = diag(3) * 0.4 + 0.05, alpha = rep(1, 3))),
    overlapping_sbm = overlapping_sbm(n = n, B = diag(3) * 0.4 + 0.05)
  )
}

test_that("step 1: closed forms equal explicit sums over a dense P", {
  for (name in names(models <- reference_models())) {
    model <- models[[name]]
    ref <- dense_reference(model)

    expect_equal(expected_edges(model), ref$edges, tolerance = 1e-10, label = name)
    expect_equal(expected_degree(model), ref$degree, tolerance = 1e-10, label = name)
    expect_equal(expected_degrees(model), unname(ref$degrees), tolerance = 1e-10, label = name)
    expect_equal(expected_density(model), ref$density, tolerance = 1e-10, label = name)

    # and the pieces are consistent with each other
    n <- model$n
    expect_equal(
      expected_density(model) * choose(n, 2),
      ref$off_diagonal_edges,
      tolerance = 1e-10,
      label = name
    )
  }
})

test_that("step 2: every sampled graph satisfies the edge counting identities", {
  models <- reference_models()

  for (name in names(models)) {
    for (loops in c(TRUE, FALSE)) {
      for (poisson in c(TRUE, FALSE)) {
        model <- models[[name]]
        model$allow_self_loops <- loops
        model$poisson_edges <- poisson
        label <- paste(name, "loops", loops, "poisson", poisson)

        for (seed in 1:5) {
          # the same seed gives the same graph in all three representations
          set.seed(seed)
          el <- sample_edgelist(model)
          set.seed(seed)
          A <- sample_sparse(model)
          set.seed(seed)
          g <- sample_igraph(model)

          tr <- sum(Matrix::diag(A))
          loops_in_el <- sum(el$from == el$to)

          # one edgelist row per undirected edge, stored with from <= to
          expect_true(all(el$from <= el$to), label = label)

          # edges between distinct nodes appear twice in A, loops once
          expect_equal(nrow(el), (sum(A) + tr) / 2, label = label)
          expect_equal(nrow(el) - loops_in_el, (sum(A) - tr) / 2, label = label)
          expect_equal(loops_in_el, tr, label = label)

          # igraph agrees on the number of edges and self-loops
          expect_equal(igraph::ecount(g), nrow(el), label = label)
          expect_equal(sum(igraph::which_loop(g)), loops_in_el, label = label)

          if (!loops) {
            expect_equal(tr, 0, label = label)
          }
        }
      }
    }
  }
})

test_that("step 3: Poisson edge totals are exactly Poisson(expected_edges())", {
  # Monte Carlo checks take ~20s in total
  skip_on_cran()

  set.seed(27)
  reps <- 400
  models <- reference_models()

  for (name in names(models)) {
    model <- models[[name]]
    mu_edges <- expected_edges(model)
    mu_off <- expected_density(model) * choose(model$n, 2)

    counts <- replicate(reps, {
      el <- sample_edgelist(model)
      c(all = nrow(el), off = sum(el$from != el$to))
    })

    z_all <- poisson_z(sum(counts["all", ]), mu_edges, reps)
    z_off <- poisson_z(sum(counts["off", ]), mu_off, reps)

    # |z| > 4.5 has probability about 7e-6 under the null, per check
    expect_lt(abs(z_all), 4.5, label = paste(name, "all edges z"))
    expect_lt(abs(z_off), 4.5, label = paste(name, "off-diagonal edges z"))

    # Poisson: variance equals mean. check the dispersion of the total
    # edge count too (ratio of sample variance to mean is ~ 1 +/- 0.07)
    dispersion <- stats::var(counts["all", ]) / mean(counts["all", ])
    expect_equal(dispersion, 1, tolerance = 0.25, label = paste(name, "dispersion"))
  }
})

test_that("without self-loops, edges between distinct nodes are still Poisson", {
  # Monte Carlo checks take ~20s in total
  skip_on_cran()

  set.seed(27)
  reps <- 400
  models <- reference_models()

  for (name in c("dense_asymmetric_S", "dcsbm", "chung_lu")) {
    model <- models[[name]]
    model$allow_self_loops <- FALSE
    mu_off <- expected_density(model) * choose(model$n, 2)

    total <- sum(replicate(reps, nrow(sample_edgelist(model))))
    expect_lt(abs(poisson_z(total, mu_off, reps)), 4.5, label = name)
  }
})

test_that("expected_density = d gives undirected graphs of density d", {
  # Monte Carlo checks take ~20s in total
  skip_on_cran()

  set.seed(27)
  n <- 400
  d <- 0.02
  reps <- 100
  pairs <- choose(n, 2)

  builders <- list(
    undirected_factor_model = function() {
      undirected_factor_model(matrix(rexp(n * 2), n), diag(2) + 0.1, expected_density = d)
    },
    erdos_renyi = function() erdos_renyi(n = n, expected_density = d),
    sbm = function() sbm(n = n, B = diag(3) * 0.5 + 0.1, expected_density = d),
    dcsbm = function() dcsbm(theta = rexp(n) + 0.5, B = diag(3) * 0.5 + 0.1, expected_density = d),
    planted_partition = function() {
      planted_partition(n = n, k = 3, within_block = 0.5, between_block = 0.1, expected_density = d)
    },
    chung_lu = function() chung_lu(theta = rexp(n) + 1, expected_density = d),
    mmsbm = function() {
      suppressMessages(mmsbm(theta = rexp(n) + 0.5, B = diag(3) * 0.5 + 0.1, alpha = rep(1, 3), expected_density = d))
    },
    overlapping_sbm = function() overlapping_sbm(n = n, B = diag(3) * 0.5 + 0.1, expected_density = d)
  )

  for (name in names(builders)) {
    model <- builders[[name]]()
    expect_equal(expected_density(model), d, label = name)

    off <- sum(replicate(reps, {
      el <- sample_edgelist(model)
      sum(el$from != el$to)
    }))

    # exact Poisson test of the realized density against d
    expect_lt(abs(poisson_z(off, d * pairs, reps)), 4.5, label = name)

    # and directly, the realized density is d (not d / 2)
    expect_equal(off / reps / pairs, d, tolerance = 0.05, label = name)
  }
})

test_that("Bernoulli SBMs realize the requested density exactly", {
  # Monte Carlo checks take ~20s in total
  skip_on_cran()

  # for blockmodels, poisson_edges = FALSE gives P(edge {i, j}) = B[z_i, z_j]
  # exactly, after rescaling B to the requested density
  set.seed(27)
  n <- 300
  d <- 0.05
  reps <- 100
  pairs <- choose(n, 2)

  model <- sbm(
    n = n,
    B = matrix(c(0.6, 0.1, 0.1, 0.4), 2),
    block_sizes = c(120, 180),
    poisson_edges = FALSE,
    allow_self_loops = FALSE,
    expected_density = d
  )

  off <- sum(replicate(reps, nrow(sample_edgelist(model))))

  # a sum of independent Bernoullis has variance at most its mean, so the
  # Poisson z-score is conservative here
  expect_lt(abs(poisson_z(off, d * pairs, reps)), 4.5)
  expect_equal(off / reps / pairs, d, tolerance = 0.03)
})

test_that("directed expectations are unchanged and exactly Poisson", {
  # Monte Carlo checks take ~20s in total
  skip_on_cran()

  set.seed(27)
  n <- 50
  d <- 40
  reps <- 400

  model <- directed_factor_model(
    X = matrix(rexp(n * 2), n),
    S = matrix(runif(6), 2, 3),
    Y = matrix(rexp(d * 3), d),
    expected_density = 0.05
  )

  P <- as.matrix(model$X %*% model$S %*% Matrix::t(model$Y))
  expect_equal(expected_edges(model), sum(P), tolerance = 1e-10)
  expect_equal(expected_density(model), mean(P), tolerance = 1e-10)
  expect_equal(expected_out_degree(model), mean(rowSums(P)), tolerance = 1e-10)
  expect_equal(expected_in_degree(model), mean(colSums(P)), tolerance = 1e-10)

  total <- sum(replicate(reps, nrow(sample_edgelist(model))))
  expect_lt(abs(poisson_z(total, expected_edges(model), reps)), 4.5)
})

test_that("erdos_renyi() requires p or a scaling argument", {
  expect_error(erdos_renyi(n = 10), "Must specify either")
  expect_silent(erdos_renyi(n = 10, expected_degree = 2))
  expect_silent(erdos_renyi(n = 10, expected_density = 0.1))

  expect_error(directed_erdos_renyi(n = 10), "Must specify either")
  expect_silent(directed_erdos_renyi(n = 10, expected_in_degree = 2))
  expect_silent(directed_erdos_renyi(n = 10, expected_density = 0.1))
})
