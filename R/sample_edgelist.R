#' Sample a random edgelist from a random dot product graph
#'
#' There are two steps to using the `fastRG` package. First,
#' you must parameterize a random dot product graph by
#' sampling the latent factors. Use functions such as
#' [dcsbm()], [sbm()], etc, to perform this specification.
#' Then, use `sample_*()` functions to generate a random graph
#' in your preferred format.
#'
#' @param factor_model A [directed_factor_model()] or
#'   [undirected_factor_model()].
#'
#' @param ... Ignored. Do not use.
#'
#' @param multigraph Logical scalar, only relevant for Poisson graphs
#'   (`poisson_edges = TRUE`), where a pair of nodes can be connected by
#'   more than one edge. When `TRUE` (the default), such repeated edges are
#'   represented as a multigraph, that is, as repeated edges (or repeated
#'   rows of an edgelist) without any edge weights. When `FALSE`, repeated
#'   edges between the same pair of nodes are collapsed into a single edge
#'   with an integer `weight` attribute that counts the number of times
#'   the edge was sampled. This weight is the corresponding entry of the
#'   adjacency matrix returned by [sample_sparse()]. In undirected graphs,
#'   the edges `(i, j)` and `(j, i)` are the same edge and are collapsed
#'   together. When `multigraph = FALSE` the `weight` attribute is always
#'   present, even for Bernoulli graphs (`poisson_edges = FALSE`), where
#'   every weight is `1`.
#'
#' @return A single realization of a random Poisson (or Bernoulli)
#'   Dot Product Graph, represented as a [tibble::tibble()] with two
#'   integer columns, `from` and `to`. When `multigraph = FALSE`, there
#'   is an additional integer column `weight` and each `(from, to)` pair
#'   appears at most once.
#'
#'   **NOTE**: Indices for isolated nodes will not appear in the edgelist!
#'   This can lead to issues if you construct network objects from the
#'   edgelist directly.
#'
#'   In the undirected case, `from` and `to` do not encode
#'   information about edge direction, but we will always have
#'   `from <= to` for convenience of edge identification.
#'
#'   To avoid handling such considerations yourself, we recommend using
#'   [sample_sparse()], [sample_igraph()], and [sample_tidygraph()]
#'   over [sample_edgelist()].
#'
#' @export
#' @family samplers
#'
#' @details This function implements the `fastRG` algorithm as
#'   described in Rohe et al (2017). Please see the paper
#'   (which is short and open access!!) for details.
#'
#' @references Rohe, Karl, Jun Tao, Xintian Han, and Norbert Binkiewicz. 2017.
#'    "A Note on Quickly Sampling a Sparse Matrix with Low Rank Expectation."
#'    Journal of Machine Learning Research; 19(77):1-13, 2018.
#'    <https://www.jmlr.org/papers/v19/17-128.html>
#'
#' @examples
#'
#' library(igraph)
#' library(tidygraph)
#'
#' set.seed(27)
#'
#' ##### undirected examples ----------------------------
#'
#' n <- 100
#' k <- 5
#'
#' X <- matrix(rpois(n = n * k, 1), nrow = n)
#' S <- matrix(runif(n = k * k, 0, .1), nrow = k)
#'
#' # S will be symmetrized internal here, or left unchanged if
#' # it is already symmetric
#'
#' ufm <- undirected_factor_model(
#'   X, S,
#'   expected_density = 0.1
#' )
#'
#' ufm
#'
#' ### sampling graphs as edgelists ----------------------
#'
#' edgelist <- sample_edgelist(ufm)
#' edgelist
#'
#' ### sampling graphs as sparse matrices ----------------
#'
#' A <- sample_sparse(ufm)
#'
#' inherits(A, "dsCMatrix")
#' isSymmetric(A)
#' dim(A)
#'
#' B <- sample_sparse(ufm)
#'
#' inherits(B, "dsCMatrix")
#' isSymmetric(B)
#' dim(B)
#'
#' ### sampling graphs as igraph graphs ------------------
#'
#' sample_igraph(ufm)
#'
#' ### sampling graphs as tidygraph graphs ---------------
#'
#' sample_tidygraph(ufm)
#'
#' ### weighted graphs instead of multigraphs ------------
#'
#' # collapse repeated edges into a single edge with an integer weight
#'
#' sample_edgelist(ufm, multigraph = FALSE)
#' sample_igraph(ufm, multigraph = FALSE)
#'
#' ##### directed examples ----------------------------
#'
#' n2 <- 100
#'
#' k1 <- 5
#' k2 <- 3
#'
#' d <- 50
#'
#' X <- matrix(rpois(n = n2 * k1, 1), nrow = n2)
#' S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
#' Y <- matrix(rexp(n = k2 * d, 1), nrow = d)
#'
#' fm <- directed_factor_model(X, S, Y, expected_in_degree = 2)
#' fm
#'
#' ### sampling graphs as edgelists ----------------------
#'
#' edgelist2 <- sample_edgelist(fm)
#' edgelist2
#'
#' ### sampling graphs as sparse matrices ----------------
#'
#' A2 <- sample_sparse(fm)
#'
#' inherits(A2, "dgCMatrix")
#' isSymmetric(A2)
#' dim(A2)
#'
#' B2 <- sample_sparse(fm)
#'
#' inherits(B2, "dgCMatrix")
#' isSymmetric(B2)
#' dim(B2)
#'
#' ### sampling graphs as igraph graphs ------------------
#'
#' # since the number of rows and the number of columns
#' # in `fm` differ, we will get a bipartite igraph here
#'
#' # creating the bipartite igraph is slow relative to other
#' # sampling -- if this is a blocker for
#' # you please open an issue and we can investigate speedups
#'
#' dig <- sample_igraph(fm)
#' is_bipartite(dig)
#'
#' ### sampling graphs as tidygraph graphs ---------------
#'
#' sample_tidygraph(fm)
#'
sample_edgelist <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  rlang::check_dots_unnamed()
  check_multigraph(multigraph)
  UseMethod("sample_edgelist")
}

#' @rdname sample_edgelist
#' @export
sample_edgelist.undirected_factor_model <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  X <- factor_model$X
  n <- factor_model$n

  # see #43, move the scaling factor here so it's temporary during sampling
  # E(A) = U S U' holds for undirected factor models
  S <- factor_model$S / 2

  # deduplication and self-loop removal must happen after the extra
  # self-loops are added below, so defer them
  # as_tibble() because the zero edge case returns a matrix
  edgelist <- tibble::as_tibble(sample_edgelist(X, S, X, FALSE, TRUE, TRUE))

  if (factor_model$allow_self_loops) {
    # off the diagonal, {i, j} gets sampled as both (i, j) and (j, i), so
    # halving S gives the right rate. (i, i) only gets sampled once, so the
    # diagonal is left at half rate. top it up with independent Poisson
    # self-loops at the missing rate diag(X (S / 2) X'), so that
    # E(A_ii) = (X S X')_ii as well. see #46
    diag_rates <- as.numeric(rowSums((X %*% S) * X))
    loops <- rep(seq_len(n), times = rpois(n, diag_rates))

    if (length(loops) > 0) {
      edgelist <- dplyr::bind_rows(edgelist, tibble(from = loops, to = loops))
    }
  }

  if (!factor_model$poisson_edges) {
    edgelist <- dplyr::distinct(edgelist)
  }

  if (!factor_model$allow_self_loops) {
    edgelist <- dplyr::filter(edgelist, to != from)
  }

  # from <= to for every edge here (including the self-loops added above),
  # so (i, j) and (j, i) are already the same row and collapsing by
  # (from, to) counts undirected edges, matching sample_sparse()
  if (!multigraph) {
    edgelist <- collapse_multiedges(edgelist)
  }

  edgelist
}

#' @rdname sample_edgelist
#' @export
sample_edgelist.directed_factor_model <- function(
  factor_model,
  ...,
  multigraph = TRUE
) {
  X <- factor_model$X
  S <- factor_model$S
  Y <- factor_model$Y

  sample_edgelist(
    X,
    S,
    Y,
    TRUE,
    factor_model$poisson_edges,
    factor_model$allow_self_loops,
    multigraph = multigraph
  )
}

#' Low level interface to sample RPDG edgelists
#'
#' **This is a brakes-off, no safety checks interface.**
#' We strongly recommend that you do not call
#' `sample_edgelist.matrix()` unless you know what you are doing,
#' and even then, we still do not recommend it, as you will
#' bypass all typical input validation.
#' **extremely loud coughing** All those who bypass input
#' validation suffer foolishly at their own hand.
#' **extremely loud coughing**
#'
#' @param factor_model An `n` by `k1` [matrix()] or [Matrix::Matrix()]
#'   of latent node positions encoding incoming edge community membership.
#'   The `X` matrix in Rohe et al (2017). Naming differs only for
#'   consistency with the S3 generic.
#'
#' @param S A `k1` by `k2` mixing [matrix()] or [Matrix::Matrix()]. In
#'   the undirect case this is assumed to be symmetric but **we do not
#'   check that this is the case**.
#'
#' @param Y A `d` by `k2` [matrix()] or [Matrix::Matrix()] of latent
#'   node positions encoding outgoing edge community membership.
#'
#' @param directed Logical indicating whether or not the graph should be
#'   directed. When `directed = FALSE`, symmetrizes `S` internally.
#'   `Y = X` together with a symmetric `S` implies a symmetric
#'   expectation (although not necessarily an undirected graph).
#'   When `directed = FALSE`, samples a directed graph with
#'   symmetric expectation, and then adds edges until symmetry
#'   is achieved.
#'
#' @param poisson_edges Whether or not to remove duplicate edges
#'   after sampling. See Section 2.3 of Rohe et al. (2017)
#'   for some additional details. Defaults to `TRUE`.
#'
#' @param allow_self_loops Logical indicating whether or not
#'   nodes should be allowed to form edges with themselves.
#'   Defaults to `TRUE`. When `FALSE`, sampling proceeds allowing
#'   self-loops, and these are then removed after the fact.
#'   Ignored for rectangular graphs (`nrow(X) != nrow(Y)`), whose rows
#'   and columns are different sets of nodes, so that there are no
#'   self-loops.
#'
#' @param ... Ignored, for generic consistency only.
#'
#' @inheritParams sample_edgelist
#' @inherit sample_edgelist return details references
#'
#' @export
#' @importFrom stats rpois rmultinom
#' @family samplers
#'
#' @examples
#'
#' set.seed(46)
#'
#' n <- 10000
#' d <- 1000
#'
#' k1 <- 5
#' k2 <- 3
#'
#' X <- matrix(rpois(n = n * k1, 1), nrow = n)
#' S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1)
#' Y <- matrix(rpois(n = d * k2, 1), nrow = d)
#'
#' sample_edgelist(X, S, Y, TRUE, TRUE, TRUE)
#'
sample_edgelist.matrix <- function(
  factor_model,
  S,
  Y,
  directed,
  poisson_edges,
  allow_self_loops,
  ...,
  multigraph = TRUE
) {
  X <- factor_model

  stopifnot(is.logical(directed))
  stopifnot(is.logical(poisson_edges))
  stopifnot(is.logical(allow_self_loops))
  check_multigraph(multigraph)

  n <- nrow(X)
  d <- nrow(Y)

  k1 <- ncol(X)
  k2 <- ncol(Y)

  Cx <- Diagonal(n = k1, x = colSums(X))
  Cy <- Diagonal(n = k2, x = colSums(Y))

  # passed to rmultinom, so Matrix objects will break things
  S_tilde <- as.matrix(Cx %*% S %*% Cy)
  expected_edges <- sum(S_tilde)

  m <- rpois(n = 1, lambda = expected_edges)

  if (m == 0) {
    # match the type of a non-empty edgelist: a tibble of integers
    edgelist <- tibble(from = integer(0), to = integer(0))

    if (!multigraph) {
      return(collapse_multiedges(edgelist))
    }

    return(edgelist)
  }

  # varpi in section 2.4 of Rohe et al (2017), the number of
  # edges between all pairs of blocks
  block_sizes <- matrix(
    rmultinom(n = 1, size = m, prob = S_tilde),
    nrow = k1,
    ncol = k2
  )

  # allocate space for the edgelist

  from <- integer(m)
  to_tmp <- integer(m)
  to <- integer(m)

  # intuition: each column works just like a block in a stochastic
  # block model. sample one column of X (variously denoted by U in
  # the paper) at a time

  u_block_start <- 1
  u_block_sizes <- rowSums(block_sizes)

  for (u in 1:k1) {
    if (u_block_sizes[u] > 0) {
      indices <- u_block_start:(u_block_start + u_block_sizes[u] - 1)

      # the prob argument should be \tilde X from the paper, but \tilde X
      # is just the l1 normalized version of X[, u], and sample()
      # will automatically normalize for us

      from[indices] <- sample(
        n,
        size = u_block_sizes[u],
        replace = TRUE,
        prob = X[, u]
      )

      u_block_start <- u_block_start + u_block_sizes[u]
    }
  }

  v_block_start <- 1
  v_block_sizes <- colSums(block_sizes)

  for (v in 1:k2) {
    if (v_block_sizes[v] > 0) {
      indices <- v_block_start:(v_block_start + v_block_sizes[v] - 1)

      # note same lack of \tilde Y as in the X/U case

      to_tmp[indices] <- sample(
        d,
        size = v_block_sizes[v],
        replace = TRUE,
        prob = Y[, v]
      )

      v_block_start <- v_block_start + v_block_sizes[v]
    }
  }

  # if the model is undirected, i think to and to_tmp are sufficient since
  # the U and V blocks should line up since block memberships match?

  # put to_tmp in the correct order, to match up with from

  u_block_start <- 1

  # i believe the move from the commented out line of code to the new
  # version should fix #13, but if something goes horribly wrong, revert this,
  # it's fine if #13 remains buggy
  # v_block_start <- c(1, cumsum(v_block_sizes))
  v_block_start <- cumsum(c(1, v_block_sizes))

  for (u in 1:k1) {
    for (v in 1:k2) {
      if (block_sizes[u, v] > 0) {
        to_index <- u_block_start:(u_block_start + block_sizes[u, v] - 1)
        tmp_index <- v_block_start[v]:(v_block_start[v] + block_sizes[u, v] - 1)

        to[to_index] <- to_tmp[tmp_index]

        v_block_start[v] <- v_block_start[v] + block_sizes[u, v]
        u_block_start <- u_block_start + block_sizes[u, v]
      }
    }
  }

  if (directed) {
    edgelist <- tibble(from = from, to = to)
  } else {
    # in the undirected case, sort the indices so that the *directed*
    # representations lives all in the same triangle (upper or lower i
    # didn't work it out)

    # *do not* move these into the tibble call otherwise you'll run
    # into a nasty NSE scoping issue that is super hard to detect
    tibble_from <- pmin(from, to)
    tibble_to <- pmax(from, to)

    edgelist <- tibble(
      from = tibble_from,
      to = tibble_to
    )
  }

  if (!poisson_edges) {
    # need to deduplicate edgelist. the number of times a given
    # (to, from) pair appears in the edgelist is the weight of
    # that edge (i.e. we're really working with a multigraph)

    edgelist <- dplyr::distinct(edgelist)
  }

  # in rectangular graphs (n != d) rows and columns are different sets of
  # nodes, so `from == to` is not a self-loop and must be kept
  if (!allow_self_loops && n == d) {
    edgelist <- dplyr::filter(edgelist, to != from)
  }

  if (!multigraph) {
    edgelist <- collapse_multiedges(edgelist)
  }

  edgelist
}

#' @rdname sample_edgelist.matrix
#' @export
sample_edgelist.Matrix <- sample_edgelist.matrix
