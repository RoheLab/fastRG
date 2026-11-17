# Create an undirected overlapping stochastic blockmodel object

To specify an overlapping stochastic blockmodel, you must specify the
number of nodes (via `n`), the mixing matrix (via `k` or `B`), and the
block membership probabilities (optional, via `pi`). In an overlapping
stochastic blockmodel, each node can belong to any number of blocks:
zero, one, or several. We provide defaults for most of these options to
enable rapid exploration, or you can invest the effort for more control
over the model parameters. We **strongly recommend** setting the
`expected_degree` or `expected_density` argument to avoid large memory
allocations associated with sampling large, dense graphs.

## Usage

``` r
overlapping_sbm(
  n,
  k = NULL,
  B = NULL,
  ...,
  pi = rep(1/k, k),
  sort_nodes = TRUE,
  force_pure = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
)
```

## Arguments

- n:

  The number of nodes in the overlapping SBM. Must be at least `k` when
  `force_pure = TRUE`.

- k:

  (mixing matrix) The number of blocks in the blockmodel. Use when you
  don't want to specify the mixing-matrix by hand. When `k` is
  specified, `B` is set to a diagonally dominant matrix with value `0.8`
  along the diagonal and `0.1 / (k - 1)` on the off-diagonal. This is
  subject to change, and may not be reproducible. `k` defaults to
  `NULL`. You must specify either `k` or `B`, but not both.

- B:

  (mixing matrix) A `k` by `k` matrix of baseline connection intensities
  between blocks. Before any rescaling via `expected_degree` or
  `expected_density`, `B[a, b]` is the expected number of edges
  contributed to the pair of nodes \\(i, j)\\ by node \\i\\'s membership
  in block `a` and node \\j\\'s membership in block `b`. These
  contributions add up over all pairs of blocks that `i` and `j` belong
  to (see "Edge formulation" below). All elements of `B` must be
  contained in \\\[0, 1\]\\. `matrix` and `Matrix` objects are both
  acceptable. If `B` is not symmetric, it will be symmetrized via the
  update `B := (B + t(B)) / 2`, and the symmetrized `B` must be
  *invertible*. Defaults to `NULL`. You must specify either `k` or `B`,
  but not both.

- ...:

  Arguments passed on to
  [`undirected_factor_model`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)

  `expected_degree`

  :   If specified, the desired expected degree of the graph. Specifying
      `expected_degree` simply rescales `S` to achieve this. Defaults to
      `NULL`. Do not specify both `expected_degree` and
      `expected_density` at the same time.

  `expected_density`

  :   If specified, the desired expected density of the graph, that is,
      the expected number of edges between distinct nodes divided by
      `choose(n, 2)`. See
      [`expected_density()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md).
      Specifying `expected_density` simply rescales `S` to achieve this.
      Defaults to `NULL`. Do not specify both `expected_degree` and
      `expected_density` at the same time.

- pi:

  (block membership probabilities) A vector of length `k` whose `j`th
  element is the probability that a node belongs to block `j`. All
  elements must be contained in \\\[0, 1\]\\. Each node joins each block
  independently, so `pi` is **not** normalized and does not need to sum
  to one: `sum(pi)` is the expected number of blocks a node joins before
  nodes in no block are reassigned (ignoring `force_pure`). A node that
  joins no block is instead assigned to exactly one block, chosen with
  probability proportional to `pi`, so every node belongs to at least
  one block and at least one element of `pi` must be positive. Defaults
  to `rep(1 / k, k)`, in which case a fraction of roughly \\(1 - 1/k)^k
  \approx e^{-1}\\ of nodes initially join no block and are reassigned
  to a single block chosen uniformly at random. Regardless of
  `sort_nodes`, `pi` is sorted in increasing order, and the rows and
  columns of `B` are permuted to match the new order of `pi`.

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block. Useful for plotting. Defaults to `TRUE`. When
  `TRUE`, the rows of `Z` are sorted in decreasing lexicographic order,
  so that members of the first block come first, and within the members
  of the first block, members of the second block come first, and so on.

- force_pure:

  Logical indicating whether or not to force presence of "pure nodes"
  (nodes that belong only to a single block) for the sake of
  identifiability. To include pure nodes, block membership sampling
  first proceeds as per usual. Then, after it is complete, `k` distinct
  nodes are chosen uniformly at random, and the memberships of the `j`th
  chosen node are overwritten so that it belongs to block `j` and to no
  other block. Requires `n >= k`. Defaults to `TRUE`.

- poisson_edges:

  Logical indicating whether or not multiple edges are allowed to form
  between a pair of nodes. Defaults to `TRUE`. When `FALSE`, sampling
  proceeds as usual, and duplicate edges are removed afterwards.
  Further, when `FALSE`, we assume that `S` specifies a desired
  Bernoulli edge probability, and back-transform `S` to a Poisson
  parameterization before sampling (see "Edge formulation" below), so
  all elements of `S` (after symmetrizing and any rescaling) must be at
  most `1`.

- allow_self_loops:

  Logical indicating whether or not nodes should be allowed to form
  edges with themselves. Defaults to `TRUE`. When `FALSE`, sampling
  proceeds allowing self-loops, and these are then removed after the
  fact.

## Value

An `undirected_overlapping_sbm` S3 object, a subclass of the
[`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)
with the following additional fields:

- `Z`: The block memberships of each node, as an `n` by `k` binary
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html).
  `Z[i, j]` is one if node `i` belongs to block `j` and zero otherwise.
  Rows can contain any number of ones, including none. Columns are
  ordered to match `pi` and `B`. For this model, `Z` is identical to the
  latent positions `X`.

- `pi`: Block membership probabilities for each block, sorted in
  increasing order.

- `B`: The mixing matrix, as specified by the user (or the default when
  only `k` is given), with rows and columns permuted to match the order
  of `pi`. Note that `B` is *not* symmetrized or rescaled: the
  symmetrized and rescaled matrix used for sampling is `S`.

- `sorted`: Logical indicating whether nodes are arranged by block
  membership (see `sort_nodes`).

## Generative Model

`overlapping_sbm()` implements the stochastic blockmodel with overlaps
(SBMO) of Kaufmann, Bonald and Lelarge (2016), with Poisson rather than
Bernoulli edges. There are two levels of randomness in an overlapping
stochastic blockmodel. First, for each node, we independently determine
if that node is a member of each block. This is handled by
`overlapping_sbm()`. Then, given these block memberships, we randomly
sample edges between nodes. This second operation is handled by
[`sample_edgelist()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md)
and
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md),
depending on your desired graph representation.

### Block memberships

Let \\Z \in \\0, 1\\^{n \times k}\\ be the binary block membership
matrix, and let \\z_i\\ be its \\i\\th row, so that \\Z\_{ia} = 1\\ when
node \\i\\ belongs to block \\a\\. Each entry is sampled independently
as

\$\$ Z\_{ia} \sim \text{Bernoulli}(\pi_a), \qquad i = 1, \dots, n, \quad
a = 1, \dots, k. \$\$

Unlike in
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
where each node belongs to exactly one block, or in
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
where each node splits one unit of membership among blocks, here a node
can be a full member of several blocks at once.

Every node must belong to at least one block. After the Bernoulli draws
above, any node with \\Z\_{ia} = 0\\ for every block \\a\\ is assigned
to exactly one block, chosen at random with probability \\\pi_a / \sum_b
\pi_b\\. As a result, the entries of \\Z\\ are not exactly independent,
and nodes that belong to a single block are somewhat more common than
independent Bernoulli draws alone would give. When `force_pure = TRUE`,
`k` randomly chosen nodes are then made pure members of one block each
(see "Identifiability").

### Degree heterogeneity

The overlapping SBM has no degree heterogeneity parameters
(equivalently, \\\theta_i = 1\\ for every node). Nodes that belong to
more blocks tend to have higher expected degree, but all members of the
same set of blocks have the same expected degree. See
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md) and
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md) for
models with degree heterogeneity.

### Edge formulation

Once we know the block memberships, we need one more ingredient, which
is the baseline intensity of connections between nodes in block `a` and
block `b`. This is given by the symmetric \\k \times k\\ matrix \\B\\.
Then each edge \\A\_{ij}\\ is Poisson distributed with parameter

\$\$ \lambda\_{ij} = z_i^T B z_j = \sum\_{a = 1}^k \sum\_{b = 1}^k
Z\_{ia} Z\_{jb} B\_{ab}. \$\$

That is, every pair of blocks \\(a, b)\\ with \\i\\ in block \\a\\ and
\\j\\ in block \\b\\ contributes \\B\_{ab}\\ to the expected number of
edges between \\i\\ and \\j\\, so two nodes that share more blocks
connect at a higher rate. When `expected_degree` or `expected_density`
is specified, \\B\\ is multiplied by a constant to achieve the desired
value, so only the relative magnitudes of the entries of \\B\\ matter in
that case.

In the factor model parameterization, the latent positions are \\X = Z\\
and the mixing matrix \\S\\ is the symmetrized and rescaled \\B\\, so
that \\\mathbb{E}\[A \mid Z\] = Z S Z^T\\.

When `poisson_edges = FALSE`, the symmetrized and rescaled \\B\\ is
instead treated as a matrix of Bernoulli probabilities, and replaced by
\\-\log(1 - B)\\ (elementwise) before sampling. Nodes \\i\\ and \\j\\
are then connected with probability

\$\$ 1 - \prod\_{a = 1}^k \prod\_{b = 1}^k (1 - B\_{ab})^{Z\_{ia}
Z\_{jb}}, \$\$

that is, each pair of shared blocks independently gets a chance to
connect \\i\\ and \\j\\. This is the "noisy-OR" model discussed in
Section 3 of Kaufmann, Bonald and Lelarge (2016), and for sparse graphs
it is close to the SBMO with Bernoulli edges.

### Identifiability

Kaufmann, Bonald and Lelarge (2016, Theorem 1) show that the parameters
\\k\\, \\B\\ and \\Z\\ of an SBMO are identifiable, up to a permutation
of the block labels, when:

1.  \\B\\ is invertible, and

2.  there is at least one "pure node" in each block that belongs to no
    other blocks.

`overlapping_sbm()` checks that the symmetrized and rescaled \\B\\ (i.e.
`S`) is invertible, and `force_pure = TRUE` guarantees the second
condition. Kaufmann, Bonald and Lelarge (2016) further assume that every
node belongs to at least one block, which `overlapping_sbm()` guarantees
by reassigning nodes that join no block (see "Block memberships" above).

### Related models

Latouche, Birmelé and Ambroise (2011) introduced the overlapping
stochastic block model (OSBM), which also uses binary membership vectors
but passes them through a logistic link, so that block contributions add
up on the logit scale rather than on the scale of the edge rate. The
OSBM is not implemented in fastRG.

The overlapping continuous community assignment model (OCCAM) of Zhang,
Levina and Zhu (2020) generalizes the SBMO by allowing non-negative,
real-valued membership vectors (normalized to have unit length) and
per-node degree heterogeneity parameters \\\theta_i\\, with
\\\mathbb{E}\[A \mid Z, \theta\] = \Theta Z B Z^T \Theta\\ for \\\Theta
= \mathrm{diag}(\theta)\\. The SBMO can be viewed as a special case of
the OCCAM (Kaufmann, Bonald and Lelarge 2016). `overlapping_sbm()` does
not include degree heterogeneity parameters.

## References

Kaufmann, Emilie, Thomas Bonald, and Marc Lelarge. 2016. "A Spectral
Algorithm with Additive Clustering for the Recovery of Overlapping
Communities in Networks." In *Algorithmic Learning Theory (ALT 2016)*,
edited by R. Ortner et al., 355–70. Lecture Notes in Computer Science,
vol. 9925. Cham: Springer International Publishing.
[doi:10.1007/978-3-319-46379-7_24](https://doi.org/10.1007/978-3-319-46379-7_24)
.

Latouche, Pierre, Etienne Birmelé, and Christophe Ambroise. 2011.
"Overlapping Stochastic Block Models with Application to the French
Political Blogosphere." *The Annals of Applied Statistics* 5 (1):
309–36. [doi:10.1214/10-AOAS382](https://doi.org/10.1214/10-AOAS382) .

Zhang, Yuan, Elizaveta Levina, and Ji Zhu. 2020. "Detecting Overlapping
Communities in Networks Using Spectral Methods." *SIAM Journal on
Mathematics of Data Science* 2 (2): 265–83.
[doi:10.1137/19M1272238](https://doi.org/10.1137/19M1272238) . Preprint:
<https://arxiv.org/abs/1412.3432>.

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md),
[`directed_planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/directed_planted_partition.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

Other undirected graphs:
[`chung_lu()`](https://rohelab.github.io/fastRG/dev/reference/chung_lu.md),
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/erdos_renyi.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

## Examples

``` r

set.seed(27)

lazy_overlapping_sbm <- overlapping_sbm(n = 1000, k = 5, expected_density = 0.01)
#> Setting `B` to a matrix with value 0.8 on the diagonal and 0.1 / (k - 1) on the off-diagonal. This parameterization may change in the future. Explicitly set `B` for reproducible results.
lazy_overlapping_sbm
#> Undirected Overlapping Stochastic Blockmodel
#> --------------------------------------------
#> 
#> Nodes (n): 1000 (sorted)
#> Blocks (k): 5
#> 
#> Traditional Overlapping SBM parameterization:
#> 
#> Block memberships (Z): 1000 x 5 [dgCMatrix] 
#> Block probabilities (pi): 5 [numeric] 
#> 
#> Block connection propensities (B): 5 x 5 [matrix] 
#> 
#> Factor model parameterization:
#> 
#> X: 1000 x 5 [dgCMatrix] 
#> S: 5 x 5 [dsyMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 5029
#> Expected degree: 10
#> Expected density: 0.01

# block memberships: a binary matrix, nodes can belong to several blocks
lazy_overlapping_sbm$Z[1:5, ]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>      block1 block2 block3 block4 block5
#> [1,]      1      1      1      1      .
#> [2,]      1      1      1      1      .
#> [3,]      1      1      1      .      1
#> [4,]      1      1      1      .      .
#> [5,]      1      1      1      .      .

# how many blocks does each node belong to? every node belongs to
# at least one block
table(Matrix::rowSums(lazy_overlapping_sbm$Z))
#> 
#>   1   2   3   4 
#> 734 186  74   6 

# sometimes you gotta let the world burn and
# sample a wildly dense graph

dense_lazy_overlapping_sbm <- overlapping_sbm(n = 500, k = 3, expected_density = 0.8)
#> Setting `B` to a matrix with value 0.8 on the diagonal and 0.1 / (k - 1) on the off-diagonal. This parameterization may change in the future. Explicitly set `B` for reproducible results.
dense_lazy_overlapping_sbm
#> Undirected Overlapping Stochastic Blockmodel
#> --------------------------------------------
#> 
#> Nodes (n): 500 (sorted)
#> Blocks (k): 3
#> 
#> Traditional Overlapping SBM parameterization:
#> 
#> Block memberships (Z): 500 x 3 [dgCMatrix] 
#> Block probabilities (pi): 3 [numeric] 
#> 
#> Block connection propensities (B): 3 x 3 [matrix] 
#> 
#> Factor model parameterization:
#> 
#> X: 500 x 3 [dgCMatrix] 
#> S: 3 x 3 [dsyMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 100645
#> Expected degree: 400.9
#> Expected density: 0.8

# explicitly setting the mixing matrix and block membership
# probabilities rather than using defaults. pi does not need to
# sum to one: here nodes belong to 1.7 blocks on average

k <- 5
n <- 1000
B <- matrix(stats::runif(k * k), nrow = k, ncol = k)

pi <- c(1, 2, 4, 1, 0.5) / 5

custom_overlapping_sbm <- overlapping_sbm(
  n = n,
  B = B,
  pi = pi,
  expected_degree = 5
)

custom_overlapping_sbm
#> Undirected Overlapping Stochastic Blockmodel
#> --------------------------------------------
#> 
#> Nodes (n): 1000 (sorted)
#> Blocks (k): 5
#> 
#> Traditional Overlapping SBM parameterization:
#> 
#> Block memberships (Z): 1000 x 5 [dgCMatrix] 
#> Block probabilities (pi): 5 [numeric] 
#> 
#> Block connection propensities (B): 5 x 5 [matrix] 
#> 
#> Factor model parameterization:
#> 
#> X: 1000 x 5 [dgCMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 2503
#> Expected degree: 5
#> Expected density: 0.005

# pi is sorted, and B is permuted to match
custom_overlapping_sbm$pi
#> [1] 0.1 0.2 0.2 0.4 0.8
colMeans(custom_overlapping_sbm$Z)
#> block1 block2 block3 block4 block5 
#>  0.099  0.215  0.211  0.394  0.843 

edgelist <- sample_edgelist(custom_overlapping_sbm)
edgelist
#> # A tibble: 2,458 × 2
#>     from    to
#>    <int> <int>
#>  1     7    82
#>  2    24    46
#>  3     6    79
#>  4    55    97
#>  5     7    11
#>  6    22    79
#>  7     4    42
#>  8    27    72
#>  9     8    10
#> 10    15    55
#> # ℹ 2,448 more rows

A <- sample_sparse(custom_overlapping_sbm)

# nodes that belong to more blocks have higher degree on average
memberships <- Matrix::rowSums(custom_overlapping_sbm$Z)
tapply(Matrix::rowSums(A), memberships, mean)
#>         1         2         3         4 
#>  2.853488  5.752525  8.724359 11.222222 

# efficient eigendecompostion that leverages low-rank structure in
# E(A) so that you don't have to form E(A) to find eigenvectors,
# as E(A) is typically dense. computation is
# handled via RSpectra

population_eigs <- eigs_sym(custom_overlapping_sbm)
population_eigs$values
#> [1]  5.94938526  0.31139205  0.10099027  0.06603763 -0.31278692
```
