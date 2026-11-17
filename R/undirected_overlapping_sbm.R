new_undirected_overlapping_sbm <- function(
  X,
  S,
  Z,
  poisson_edges = TRUE,
  allow_self_loops = TRUE,
  pi,
  B,
  sorted,
  ...,
  subclass = character()
) {
  subclass <- c(subclass, "undirected_overlapping_sbm")
  overlapping_sbm <- undirected_factor_model(X, S, ..., subclass = subclass)
  overlapping_sbm$Z <- Z
  overlapping_sbm$pi <- pi
  overlapping_sbm$B <- B
  overlapping_sbm$sorted <- sorted
  overlapping_sbm
}

validate_undirected_overlapping_sbm <- function(x) {
  values <- unclass(x)

  if (any(values$pi < 0) || any(values$pi > 1)) {
    stop(
      "All elements of `pi` must be contained in [0, 1].",
      call. = FALSE
    )
  }

  if (!is_invertible(values$S)) {
    stop(
      "`B` must be an invertible matrix.",
      call. = FALSE
    )
  }

  if (max(values$B) > 1 || min(values$B) < 0) {
    stop(
      "All elements of `B` must be contained in [0, 1].",
      call. = FALSE
    )
  }

  x
}

#' Create an undirected overlapping stochastic blockmodel object
#'
#' To specify an overlapping stochastic blockmodel, you must specify
#' the number of nodes (via `n`),
#' the mixing matrix (via `k` or `B`), and the block membership
#' probabilities (optional, via `pi`). In an overlapping stochastic
#' blockmodel, each node can belong to any number of blocks: zero, one,
#' or several. We provide defaults for most of these
#' options to enable rapid exploration, or you can invest the effort
#' for more control over the model parameters. We **strongly recommend**
#' setting the `expected_degree` or `expected_density` argument
#' to avoid large memory allocations associated with
#' sampling large, dense graphs.
#'
#' @param n The number of nodes in the overlapping SBM. Must be at least
#'   `k` when `force_pure = TRUE`.
#'
#' @param k (mixing matrix) The number of blocks in the blockmodel.
#'   Use when you don't want to specify the mixing-matrix by hand.
#'   When `k` is specified, `B` is set to a diagonally dominant matrix with
#'   value `0.8` along the diagonal and `0.1 / (k - 1)` on the
#'   off-diagonal. This is subject to change, and may not be reproducible.
#'   `k` defaults to `NULL`. You must specify either `k`
#'   or `B`, but not both.
#'
#' @param B (mixing matrix) A `k` by `k` matrix of baseline connection
#'   intensities between blocks. Before any rescaling via
#'   `expected_degree` or `expected_density`, `B[a, b]` is the expected
#'   number of edges contributed to the pair of nodes \eqn{(i, j)} by
#'   node \eqn{i}'s membership in block `a` and node \eqn{j}'s membership
#'   in block `b`. These contributions add up over all pairs of blocks
#'   that `i` and `j` belong to (see "Edge formulation" below).
#'   All elements of `B` must be contained in \eqn{[0, 1]}.
#'   `matrix` and `Matrix` objects are both
#'   acceptable. If `B` is not symmetric, it will be
#'   symmetrized via the update `B := (B + t(B)) / 2`, and the
#'   symmetrized `B` must be *invertible*. Defaults to `NULL`.
#'   You must specify either `k` or `B`, but not both.
#'
#' @param pi (block membership probabilities) A vector of length `k`
#'   whose `j`th element is the probability that a node belongs to
#'   block `j`. All elements must be contained in \eqn{[0, 1]}. Each node
#'   joins each block independently, so `pi` is **not** normalized and
#'   does not need to sum to one: `sum(pi)` is the expected number of
#'   blocks a node joins before nodes in no block are reassigned
#'   (ignoring `force_pure`). A node that joins no block is instead
#'   assigned to exactly one block, chosen with probability
#'   proportional to `pi`, so every node belongs to at least one block
#'   and at least one element of `pi` must be positive. Defaults to
#'   `rep(1 / k, k)`, in which case a fraction of roughly
#'   \eqn{(1 - 1/k)^k \approx e^{-1}}{(1 - 1/k)^k, approximately exp(-1),}
#'   of nodes initially join no block and are reassigned to a single
#'   block chosen uniformly at random. Regardless of `sort_nodes`,
#'   `pi` is sorted in increasing order, and the rows and columns
#'   of `B` are permuted to match the new order of `pi`.
#'
#' @param sort_nodes Logical indicating whether or not to sort the nodes
#'   so that they are grouped by block. Useful for plotting.
#'   Defaults to `TRUE`. When `TRUE`, the rows of `Z` are sorted in
#'   decreasing lexicographic order, so that members of the first
#'   block come first, and within the members of the first block,
#'   members of the second block come first, and so on.
#'
#' @param force_pure Logical indicating whether or not to force presence of
#'   "pure nodes" (nodes that belong only to a single block) for the sake
#'   of identifiability. To include pure nodes, block membership sampling
#'   first proceeds as per usual. Then, after it is complete, `k` distinct
#'   nodes are chosen uniformly at random, and the memberships of the
#'   `j`th chosen node are overwritten so that it belongs to block `j`
#'   and to no other block. Requires `n >= k`. Defaults to `TRUE`.
#'
#' @param poisson_edges Currently **ignored** by `overlapping_sbm()`.
#'   The argument is accepted but not passed on to
#'   [undirected_factor_model()], so the returned object always has
#'   `poisson_edges = TRUE`, `S` is not back-transformed to a Bernoulli
#'   parameterization, and sampled graphs may contain multiple edges
#'   between a pair of nodes.
#'
#' @param allow_self_loops Currently **ignored** by `overlapping_sbm()`.
#'   The argument is accepted but not passed on to
#'   [undirected_factor_model()], so the returned object always has
#'   `allow_self_loops = TRUE` and sampled graphs may contain self-loops.
#'
#' @inheritDotParams undirected_factor_model expected_degree expected_density
#'
#' @return An `undirected_overlapping_sbm` S3 object, a subclass of the
#'   [undirected_factor_model()] with the following additional
#'   fields:
#'
#'   - `Z`: The block memberships of each node, as an `n` by `k`
#'     binary [Matrix::Matrix()]. `Z[i, j]` is one if node `i` belongs
#'     to block `j` and zero otherwise. Rows can contain any number of
#'     ones, including none. Columns are ordered to match `pi` and `B`.
#'     For this model, `Z` is identical to the latent positions `X`.
#'
#'   - `pi`: Block membership probabilities for each block, sorted in
#'     increasing order.
#'
#'   - `B`: The mixing matrix, as specified by the user (or the default
#'     when only `k` is given), with rows and columns permuted to match
#'     the order of `pi`. Note that `B` is *not* symmetrized or rescaled:
#'     the symmetrized and rescaled matrix used for sampling is `S`.
#'
#'   - `sorted`: Logical indicating whether nodes are arranged by
#'     block membership (see `sort_nodes`).
#'
#' @references
#'
#' Kaufmann, Emilie, Thomas Bonald, and Marc Lelarge. 2016.
#' "A Spectral Algorithm with Additive Clustering for the Recovery of
#' Overlapping Communities in Networks." In *Algorithmic Learning
#' Theory (ALT 2016)*, edited by R. Ortner et al., 355–70.
#' Lecture Notes in Computer Science, vol. 9925.
#' Cham: Springer International Publishing.
#' \doi{10.1007/978-3-319-46379-7_24}.
#'
#' Latouche, Pierre, Etienne Birmelé, and Christophe Ambroise. 2011.
#' "Overlapping Stochastic Block Models with Application to the
#' French Political Blogosphere." *The Annals of Applied Statistics* 5
#' (1): 309–36. \doi{10.1214/10-AOAS382}.
#'
#' Zhang, Yuan, Elizaveta Levina, and Ji Zhu. 2020. "Detecting
#' Overlapping Communities in Networks Using Spectral Methods."
#' *SIAM Journal on Mathematics of Data Science* 2 (2): 265–83.
#' \doi{10.1137/19M1272238}. Preprint: <https://arxiv.org/abs/1412.3432>.
#'
#' @export
#' @family stochastic block models
#' @family undirected graphs
#'
#' @details
#'
#' # Generative Model
#'
#' `overlapping_sbm()` implements the stochastic blockmodel with
#' overlaps (SBMO) of Kaufmann, Bonald and Lelarge (2016), with Poisson
#' rather than Bernoulli edges. There are two levels of randomness in an
#' overlapping stochastic blockmodel. First, for each node, we
#' independently determine if that node is a member of each block. This is
#' handled by `overlapping_sbm()`. Then, given these block memberships,
#' we randomly sample edges between nodes. This second
#' operation is handled by [sample_edgelist()],
#' [sample_sparse()], [sample_igraph()] and
#' [sample_tidygraph()], depending on your desired
#' graph representation.
#'
#' ## Block memberships
#'
#' Let \eqn{Z \in \{0, 1\}^{n \times k}}{Z} be the binary block
#' membership matrix, and let \eqn{z_i}{z[i]} be its \eqn{i}th row, so
#' that \eqn{Z_{ia} = 1}{Z[i, a] = 1} when node \eqn{i} belongs to block
#' \eqn{a}. Each entry is sampled independently as
#'
#' \deqn{
#'   Z_{ia} \sim \text{Bernoulli}(\pi_a),
#'   \qquad i = 1, \dots, n, \quad a = 1, \dots, k.
#' }{
#'   Z[i, a] ~ Bernoulli(pi[a]), i = 1, ..., n, a = 1, ..., k.
#' }
#'
#' Unlike in [dcsbm()], where each node belongs to exactly one block,
#' or in [mmsbm()], where each node splits one unit of membership among
#' blocks, here a node can be a full member of several blocks at
#' once.
#'
#' Every node must belong to at least one block. After the Bernoulli
#' draws above, any node with \eqn{Z_{ia} = 0}{Z[i, a] = 0} for every
#' block \eqn{a} is assigned to exactly one block, chosen at random with
#' probability \eqn{\pi_a / \sum_b \pi_b}{pi[a] / sum(pi)}. As a result,
#' the entries of \eqn{Z} are not exactly independent, and nodes that
#' belong to a single block are somewhat more common than independent
#' Bernoulli draws alone would give. When `force_pure = TRUE`, `k` randomly chosen nodes
#' are then made pure members of one block each (see "Identifiability").
#'
#' ## Degree heterogeneity
#'
#' The overlapping SBM has no degree heterogeneity parameters
#' (equivalently, \eqn{\theta_i = 1} for every node). Nodes that belong
#' to more blocks tend to have higher expected degree, but all members
#' of the same set of blocks have the same expected degree. See [dcsbm()]
#' and [mmsbm()] for models with degree heterogeneity.
#'
#' ## Edge formulation
#'
#' Once we know the block memberships, we need one more
#' ingredient, which is the baseline intensity of connections
#' between nodes in block `a` and block `b`. This is given by the
#' symmetric \eqn{k \times k}{k by k} matrix \eqn{B}. Then each edge
#' \eqn{A_{ij}}{A[i, j]} is Poisson distributed with parameter
#'
#' \deqn{
#'   \lambda_{ij} = z_i^T B z_j
#'     = \sum_{a = 1}^k \sum_{b = 1}^k Z_{ia} Z_{jb} B_{ab}.
#' }{
#'   lambda[i, j] = z[i]' B z[j] = sum_a sum_b Z[i, a] Z[j, b] B[a, b].
#' }
#'
#' That is, every pair of blocks \eqn{(a, b)} with \eqn{i} in block
#' \eqn{a} and \eqn{j} in block \eqn{b} contributes \eqn{B_{ab}}{B[a, b]} to the
#' expected number of edges between \eqn{i} and \eqn{j}, so two nodes
#' that share more blocks connect at a higher rate. When
#' `expected_degree` or `expected_density` is specified, \eqn{B} is
#' multiplied by a constant to achieve the desired value, so only the
#' relative magnitudes of the entries of \eqn{B} matter in that case.
#'
#' In the factor model parameterization, the latent positions are
#' \eqn{X = Z} and the mixing matrix \eqn{S} is the symmetrized and
#' rescaled \eqn{B}, so that
#' \eqn{\mathbb{E}[A \mid Z] = Z S Z^T}{E[A | Z] = Z S Z'}.
#'
#' ## Identifiability
#'
#' Kaufmann, Bonald and Lelarge (2016, Theorem 1) show that the
#' parameters \eqn{k}, \eqn{B} and \eqn{Z} of an SBMO are identifiable,
#' up to a permutation of the block labels, when:
#'
#' 1. \eqn{B} is invertible, and
#' 2. there is at least one "pure node" in each block that belongs to no
#'   other blocks.
#'
#' `overlapping_sbm()` checks that the symmetrized and rescaled
#' \eqn{B} (i.e. `S`) is invertible, and `force_pure = TRUE` guarantees
#' the second condition. Kaufmann, Bonald and Lelarge (2016) further
#' assume that every node belongs to at least one block, which
#' `overlapping_sbm()` guarantees by reassigning nodes that join no
#' block (see "Block memberships" above).
#'
#' ## Related models
#'
#' Latouche, Birmelé and Ambroise (2011) introduced the overlapping
#' stochastic block model (OSBM), which also uses binary membership
#' vectors but passes them through a logistic link, so that block
#' contributions add up on the logit scale rather than on the scale of
#' the edge rate. The OSBM is not implemented in fastRG.
#'
#' The overlapping continuous community assignment model (OCCAM) of Zhang,
#' Levina and Zhu (2020) generalizes the SBMO by allowing non-negative,
#' real-valued membership vectors (normalized to have unit length) and
#' per-node degree heterogeneity parameters \eqn{\theta_i}, with
#' \eqn{\mathbb{E}[A \mid Z, \theta] = \Theta Z B Z^T \Theta}{E[A | Z, theta] = Theta Z B Z' Theta}
#' for \eqn{\Theta = \mathrm{diag}(\theta)}{Theta = diag(theta)}. The
#' SBMO can be viewed as a special case of the OCCAM (Kaufmann, Bonald
#' and Lelarge 2016). `overlapping_sbm()` does not include degree
#' heterogeneity parameters.
#'
#' @examples
#'
#' set.seed(27)
#'
#' lazy_overlapping_sbm <- overlapping_sbm(n = 1000, k = 5, expected_density = 0.01)
#' lazy_overlapping_sbm
#'
#' # block memberships: a binary matrix, nodes can belong to several blocks
#' lazy_overlapping_sbm$Z[1:5, ]
#'
#' # how many blocks does each node belong to? every node belongs to
#' # at least one block
#' table(Matrix::rowSums(lazy_overlapping_sbm$Z))
#'
#' # sometimes you gotta let the world burn and
#' # sample a wildly dense graph
#'
#' dense_lazy_overlapping_sbm <- overlapping_sbm(n = 500, k = 3, expected_density = 0.8)
#' dense_lazy_overlapping_sbm
#'
#' # explicitly setting the mixing matrix and block membership
#' # probabilities rather than using defaults. pi does not need to
#' # sum to one: here nodes belong to 1.7 blocks on average
#'
#' k <- 5
#' n <- 1000
#' B <- matrix(stats::runif(k * k), nrow = k, ncol = k)
#'
#' pi <- c(1, 2, 4, 1, 0.5) / 5
#'
#' custom_overlapping_sbm <- overlapping_sbm(
#'   n = n,
#'   B = B,
#'   pi = pi,
#'   expected_degree = 5
#' )
#'
#' custom_overlapping_sbm
#'
#' # pi is sorted, and B is permuted to match
#' custom_overlapping_sbm$pi
#' colMeans(custom_overlapping_sbm$Z)
#'
#' edgelist <- sample_edgelist(custom_overlapping_sbm)
#' edgelist
#'
#' A <- sample_sparse(custom_overlapping_sbm)
#'
#' # nodes that belong to more blocks have higher degree on average
#' memberships <- Matrix::rowSums(custom_overlapping_sbm$Z)
#' tapply(Matrix::rowSums(A), memberships, mean)
#'
#' # efficient eigendecompostion that leverages low-rank structure in
#' # E(A) so that you don't have to form E(A) to find eigenvectors,
#' # as E(A) is typically dense. computation is
#' # handled via RSpectra
#'
#' population_eigs <- eigs_sym(custom_overlapping_sbm)
#' population_eigs$values
#'
overlapping_sbm <- function(
  n,
  k = NULL,
  B = NULL,
  ...,
  pi = rep(1 / k, k),
  sort_nodes = TRUE,
  force_pure = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
) {
  ### mixing matrix

  if (is.null(k) && is.null(B)) {
    stop("Must specify either `k` or `B`.", call. = FALSE)
  } else if (is.null(B)) {
    if (k < 1) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    message(
      "Setting `B` to a matrix with value 0.8 on the diagonal and ",
      "0.1 / (k - 1) on the off-diagonal. This parameterization may change ",
      "in the future. Explicitly set `B` for reproducible results."
    )

    B <- matrix(0.1 / (k - 1), nrow = k, ncol = k)
    diag(B) <- 0.8
  } else if (is.null(k)) {
    if (nrow(B) != ncol(B)) {
      stop("`B` must be a square matrix.", call. = FALSE)
    }

    k <- nrow(B)
  }

  ### block membership

  if (length(pi) != nrow(B) || length(pi) != ncol(B)) {
    stop("Length of `pi` must match dimensions of `B`.", call. = FALSE)
  }

  if (any(pi < 0) || any(pi > 1)) {
    stop("All elements of `pi` must be contained in [0, 1].", call. = FALSE)
  }

  if (sum(pi) <= 0) {
    stop("At least one element of `pi` must be positive.", call. = FALSE)
  }

  # order mixing matrix by expected group size

  if (k > 1) {
    B <- B[order(pi), ]
    B <- B[, order(pi)]
  }

  pi <- sort(pi) # do not normalize pi!

  # sample block memberships

  X <- Matrix(0, nrow = n, ncol = k)
  colnames(X) <- paste0("block", 1:k)

  for (i in 1:k) {
    ind <- stats::rbinom(n, 1, pi[i])
    X[, i] <- ind
  }

  # every node must belong to at least one block. nodes that joined no
  # block above are assigned to exactly one block, chosen with
  # probability proportional to pi

  empty <- which(Matrix::rowSums(X) == 0)

  if (length(empty) > 0) {
    blocks <- sample.int(k, length(empty), replace = TRUE, prob = pi)
    X[cbind(empty, blocks)] <- 1
  }

  if (force_pure) {
    pure_indices <- sample(n, k)

    for (i in 1:k) {
      X[pure_indices[i], ] <- 0
      X[pure_indices[i], i] <- 1
    }
  }

  if (sort_nodes) {
    X <- sort_by_all_columns(X)
  }

  overlapping_sbm <- new_undirected_overlapping_sbm(
    X = X,
    S = B, # accepts B but transforms by symmetrizing and scaling internally
    B = B,
    Z = X,
    pi = pi,
    sorted = sort_nodes,
    poisson_edges = poisson_edges,
    allow_self_loops = allow_self_loops,
    ...
  )

  validate_undirected_overlapping_sbm(overlapping_sbm)
}

#' @method print undirected_overlapping_sbm
#' @export
print.undirected_overlapping_sbm <- function(x, ...) {
  cat(glue(
    "Undirected Degree-Corrected Overlapping Blockmodel\n",
    .trim = FALSE
  ))
  cat(glue(
    "-------------------------------------------------\n\n",
    .trim = FALSE
  ))

  sorted <- if (x$sorted) "sorted" else "sorted"

  cat(glue("Nodes (n): {x$n} ({sorted})\n", .trim = FALSE))
  cat(glue("Blocks (k): {x$k}\n\n", .trim = FALSE))

  cat("Traditional Overlapping SBM parameterization:\n\n")
  cat("Block memberships (Z):", dim_and_class(x$Z), "\n")
  cat("Block probabilities (pi):", dim_and_class(x$pi), "\n\n")
  cat("Block connection propensities (B):", dim_and_class(x$B), "\n\n")

  cat("Factor model parameterization:\n\n")
  cat("X:", dim_and_class(x$X), "\n")
  cat("S:", dim_and_class(x$S), "\n\n")

  cat("Poisson edges:", as.character(x$poisson_edges), "\n")
  cat("Allow self loops:", as.character(x$allow_self_loops), "\n\n")

  cat(glue("Expected edges: {round(expected_edges(x))}\n", .trim = FALSE))
  cat(glue("Expected degree: {round(expected_degree(x), 1)}\n", .trim = FALSE))
  cat(glue("Expected density: {round(expected_density(x), 5)}", .trim = FALSE))
}
