# Create an undirected overlapping degree corrected stochastic blockmodel object

To specify a overlapping stochastic blockmodel, you must specify the
number of nodes (via `n`), the mixing matrix (via `k` or `B`), and the
block probabilities (optional, via `pi`). We provide defaults for most
of these options to enable rapid exploration, or you can invest the
effort for more control over the model parameters. We **strongly
recommend** setting the `expected_degree` or `expected_density` argument
to avoid large memory allocations associated with sampling large, dense
graphs.

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

  The number of nodes in the overlapping SBM.

- k:

  (mixing matrix) The number of blocks in the blockmodel. Use when you
  don't want to specify the mixing-matrix by hand. When `k` is
  specified, `B` is set to a diagonal dominant matrix with value `0.8`
  along the diagonal and `0.1 / (k - 1)` on the off-diagonal. `k`
  defaults to `NULL`. You must specify either `k` or `B`, but not both.

- B:

  (mixing matrix) A `k` by `k` matrix of block connection probabilities.
  The probability that a node in block `i` connects to a node in
  community `j` is `Poisson(B[i, j])`. Must be an *invertible*,
  symmetric square matrix. `matrix` and `Matrix` objects are both
  acceptable. If `B` is not symmetric, it will be symmetrized via the
  update `B := B + t(B)`. Defaults to `NULL`. You must specify either
  `k` or `B`, but not both.

- ...:

  Arguments passed on to
  [`undirected_factor_model`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)

  `expected_degree`

  :   If specified, the desired expected degree of the graph. Specifying
      `expected_degree` simply rescales `S` to achieve this. Defaults to
      `NULL`. Do not specify both `expected_degree` and
      `expected_density` at the same time.

  `expected_density`

  :   If specified, the desired expected density of the graph.
      Specifying `expected_density` simply rescales `S` to achieve this.
      Defaults to `NULL`. Do not specify both `expected_degree` and
      `expected_density` at the same time.

- pi:

  (block probabilities) Probability of membership in each block.
  Membership in each block is independent under the overlapping SBM.
  Defaults to `rep(1 / k, k)`.

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block. Useful for plotting. Defaults to `TRUE`.

- force_pure:

  Logical indicating whether or not to force presence of "pure nodes"
  (nodes that belong only to a single community) for the sake of
  identifiability. To include pure nodes, block membership sampling
  first proceeds as per usual. Then, after it is complete, `k` nodes are
  chosen randomly as pure nodes, one for each block. Defaults to `TRUE`.

- poisson_edges:

  Logical indicating whether or not multiple edges are allowed to form
  between a pair of nodes. Defaults to `TRUE`. When `FALSE`, sampling
  proceeds as usual, and duplicate edges are removed afterwards.
  Further, when `FALSE`, we assume that `S` specifies a desired
  between-factor connection probability, and back-transform this `S` to
  the appropriate Poisson intensity parameter to approximate Bernoulli
  factor connection probabilities. See Section 2.3 of Rohe et al. (2017)
  for some additional details.

- allow_self_loops:

  Logical indicating whether or not nodes should be allowed to form
  edges with themselves. Defaults to `TRUE`. When `FALSE`, sampling
  proceeds allowing self-loops, and these are then removed after the
  fact.

## Value

An `undirected_overlapping_sbm` S3 object, a subclass of the
[`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)
with the following additional fields:

- `pi`: Sampling probabilities for each block.

- `sorted`: Logical indicating where nodes are arranged by block (and
  additionally by degree heterogeneity parameter) within each block.

## Generative Model

There are two levels of randomness in a degree-corrected overlapping
stochastic blockmodel. First, for each node, we independently determine
if that node is a member of each block. This is handled by
`overlapping_sbm()`. Then, given these block memberships, we randomly
sample edges between nodes. This second operation is handled by
[`sample_edgelist()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md)
and
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md),
depending depending on your desired graph representation.

### Identifiability

In order to be identifiable, an overlapping SBM must satisfy two
conditions:

1.  `B` must be invertible, and

2.  the must be at least one "pure node" in each block that belongs to
    no other blocks.

### Block memberships

Note that some nodes may not belong to any blocks.

**TODO**

### Edge formulation

Once we know the block memberships, we need one more ingredient, which
is the baseline intensity of connections between nodes in block `i` and
block `j`. Then each edge \\A\_{i,j}\\ is Poisson distributed with
parameter

**TODO**

## References

Kaufmann, Emilie, Thomas Bonald, and Marc Lelarge. "A Spectral Algorithm
with Additive Clustering for the Recovery of Overlapping Communities in
Networks," Vol. 9925. Lecture Notes in Computer Science. Cham: Springer
International Publishing, 2016.
https://doi.org/10.1007/978-3-319-46379-7.

Latouche, Pierre, Etienne Birmelé, and Christophe Ambroise. "Overlapping
Stochastic Block Models with Application to the French Political
Blogosphere." The Annals of Applied Statistics 5, no. 1 (March 2011):
309–36. https://doi.org/10.1214/10-AOAS382.

Zhang, Yuan, Elizaveta Levina, and Ji Zhu. "Detecting Overlapping
Communities in Networks Using Spectral Methods." ArXiv:1412.3432,
December 10, 2014. http://arxiv.org/abs/1412.3432.

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md),
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
#> Undirected Degree-Corrected Overlapping Blockmodel
#> -------------------------------------------------
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
#> Expected edges: 4995
#> Expected degree: 5
#> Expected density: 0.01

# sometimes you gotta let the world burn and
# sample a wildly dense graph

dense_lazy_overlapping_sbm <- overlapping_sbm(n = 500, k = 3, expected_density = 0.8)
#> Setting `B` to a matrix with value 0.8 on the diagonal and 0.1 / (k - 1) on the off-diagonal. This parameterization may change in the future. Explicitly set `B` for reproducible results.
dense_lazy_overlapping_sbm
#> Undirected Degree-Corrected Overlapping Blockmodel
#> -------------------------------------------------
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
#> Expected edges: 99800
#> Expected degree: 199.6
#> Expected density: 0.8

k <- 5
n <- 1000
B <- matrix(stats::runif(k * k), nrow = k, ncol = k)

pi <- c(1, 2, 4, 1, 1) / 5

custom_overlapping_sbm <- overlapping_sbm(
  n = 200,
  B = B,
  pi = pi,
  expected_degree = 5
)

custom_overlapping_sbm
#> Undirected Degree-Corrected Overlapping Blockmodel
#> -------------------------------------------------
#> 
#> Nodes (n): 200 (sorted)
#> Blocks (k): 5
#> 
#> Traditional Overlapping SBM parameterization:
#> 
#> Block memberships (Z): 200 x 5 [dgCMatrix] 
#> Block probabilities (pi): 5 [numeric] 
#> 
#> Block connection propensities (B): 5 x 5 [matrix] 
#> 
#> Factor model parameterization:
#> 
#> X: 200 x 5 [dgCMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 1000
#> Expected degree: 5
#> Expected density: 0.05025

edgelist <- sample_edgelist(custom_overlapping_sbm)
edgelist
#> # A tibble: 1,026 × 2
#>     from    to
#>    <int> <int>
#>  1     3    33
#>  2     2     8
#>  3     1     1
#>  4    28    39
#>  5    26    33
#>  6    10    19
#>  7    19    32
#>  8    21    24
#>  9    25    41
#> 10    23    33
#> # ℹ 1,016 more rows

# efficient eigendecompostion that leverages low-rank structure in
# E(A) so that you don't have to form E(A) to find eigenvectors,
# as E(A) is typically dense. computation is
# handled via RSpectra

population_eigs <- eigs_sym(custom_overlapping_sbm)
```
