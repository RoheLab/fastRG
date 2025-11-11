# Create an undirected degree-corrected mixed membership stochastic blockmodel object

To specify a degree-corrected mixed membership stochastic blockmodel,
you must specify the degree-heterogeneity parameters (via `n` or
`theta`), the mixing matrix (via `k` or `B`), and the relative block
propensities (optional, via `alpha`). We provide defaults for most of
these options to enable rapid exploration, or you can invest the effort
for more control over the model parameters. We **strongly recommend**
setting the `expected_degree` or `expected_density` argument to avoid
large memory allocations associated with sampling large, dense graphs.

## Usage

``` r
mmsbm(
  n = NULL,
  theta = NULL,
  k = NULL,
  B = NULL,
  ...,
  alpha = rep(1, k),
  sort_nodes = TRUE,
  force_pure = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
)
```

## Arguments

- n:

  (degree heterogeneity) The number of nodes in the blockmodel. Use when
  you don't want to specify the degree-heterogeneity parameters `theta`
  by hand. When `n` is specified, `theta` is randomly generated from a
  `LogNormal(2, 1)` distribution. This is subject to change, and may not
  be reproducible. `n` defaults to `NULL`. You must specify either `n`
  or `theta`, but not both.

- theta:

  (degree heterogeneity) A numeric vector explicitly specifying the
  degree heterogeneity parameters. This implicitly determines the number
  of nodes in the resulting graph, i.e. it will have `length(theta)`
  nodes. Must be positive. Setting to a vector of ones recovers a
  stochastic blockmodel without degree correction. Defaults to `NULL`.
  You must specify either `n` or `theta`, but not both.

- k:

  (mixing matrix) The number of blocks in the blockmodel. Use when you
  don't want to specify the mixing-matrix by hand. When `k` is
  specified, the elements of `B` are drawn randomly from a
  `Uniform(0, 1)` distribution. This is subject to change, and may not
  be reproducible. `k` defaults to `NULL`. You must specify either `k`
  or `B`, but not both.

- B:

  (mixing matrix) A `k` by `k` matrix of block connection probabilities.
  The probability that a node in block `i` connects to a node in
  community `j` is `Poisson(B[i, j])`. Must be a square matrix. `matrix`
  and `Matrix` objects are both acceptable. If `B` is not symmetric, it
  will be symmetrized via the update `B := B + t(B) / 2`. Defaults to
  `NULL`. You must specify either `k` or `B`, but not both.

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

- alpha:

  (relative block propensities) Relative block propensities, which are
  parameters of a Dirichlet distribution. All elments of `alpha` must
  thus be positive. Must match the dimensions of `B` or `k`. Defaults to
  `rep(1, k)`, or balanced membership across blocks.

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block and by `theta`. Useful for plotting. Defaults to
  `TRUE`. When `TRUE`, nodes are first sorted by block membership, and
  then by degree-correction parameters within each block. Additionally,
  `pi` is sorted in increasing order, and the columns of the `B` matrix
  are permuted to match the new order of `pi`.

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

An `undirected_mmsbm` S3 object, a subclass of the
[`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)
with the following additional fields:

- `theta`: A numeric vector of degree-heterogeneity parameters.

- `Z`: The community memberships of each node, a
  [`matrix()`](https://rdrr.io/r/base/matrix.html) with `k` columns,
  whose row sums all equal one.

- `alpha`: Community membership proportion propensities.

- `sorted`: Logical indicating where nodes are arranged by block (and
  additionally by degree heterogeneity parameter) within each block.

## Generative Model

There are two levels of randomness in a degree-corrected stochastic
blockmodel. First, we randomly choose how much each node belongs to each
block in the blockmodel. Each node is one unit of block membership to
distribute. This is handled by `mmsbm()`. Then, given these block
memberships, we randomly sample edges between nodes. This second
operation is handled by
[`sample_edgelist()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md)
and
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md),
depending depending on your desired graph representation.

### Block memberships

Let \\Z_i\\ by a vector on the `k` dimensional simplex representing the
block memberships of node \\i\\. To generate \\z_i\\ we sample from a
Dirichlet distribution with parameter vector \\\alpha\\. Block
memberships for each node are independent.

### Degree heterogeneity

In addition to block membership, the MMSBM also allows nodes to have
different propensities for edge formation. We represent this propensity
for node \\i\\ by a positive number \\\theta_i\\.

### Edge formulation

Once we know the block membership vector \\z_i, z_j\\ and the degree
heterogeneity parameters \\\theta\\, we need one more ingredient, which
is the baseline intensity of connections between nodes in block `i` and
block `j`. This is given by a \\k \times k\\ matrix \\B\\. Then each
edge \\A\_{i,j}\\ is Poisson distributed with parameter

\$\$ \lambda\_{i, j} = \theta_i \cdot z_i^T B z_j \cdot \theta_j. \$\$

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

Other undirected graphs:
[`chung_lu()`](https://rohelab.github.io/fastRG/dev/reference/chung_lu.md),
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/erdos_renyi.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

## Examples

``` r
set.seed(27)

lazy_mmsbm <- mmsbm(n = 100, k = 5, expected_density = 0.01)
#> Generating random degree heterogeneity parameters `theta` from a LogNormal(2, 1) distribution. This distribution may change in the future. Explicitly set `theta` for reproducible results.
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.
lazy_mmsbm
#> Undirected Degree-Corrected Mixed Membership Stochastic Blockmodel
#> ------------------------------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 5
#> 
#> Traditional MMSBM parameterization:
#> 
#> Block memberships portions (Z): 100 x 5 [matrix] 
#> Degree heterogeneity (theta): 100 [numeric] 
#> Block propensities (alpha): 5 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 5 [dgeMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 50
#> Expected degree: 0.5
#> Expected density: 0.01

# sometimes you gotta let the world burn and
# sample a wildly dense graph

dense_lazy_mmsbm <- mmsbm(n = 500, k = 3, expected_density = 0.8)
#> Generating random degree heterogeneity parameters `theta` from a LogNormal(2, 1) distribution. This distribution may change in the future. Explicitly set `theta` for reproducible results.
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.
dense_lazy_mmsbm
#> Undirected Degree-Corrected Mixed Membership Stochastic Blockmodel
#> ------------------------------------------------------------------
#> 
#> Nodes (n): 500 (arranged by block)
#> Blocks (k): 3
#> 
#> Traditional MMSBM parameterization:
#> 
#> Block memberships portions (Z): 500 x 3 [matrix] 
#> Degree heterogeneity (theta): 500 [numeric] 
#> Block propensities (alpha): 3 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 500 x 3 [dgeMatrix] 
#> S: 3 x 3 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 99800
#> Expected degree: 199.6
#> Expected density: 0.8

# explicitly setting the degree heterogeneity parameter,
# mixing matrix, and relative community sizes rather
# than using randomly generated defaults

k <- 5
n <- 100
B <- matrix(stats::runif(k * k), nrow = k, ncol = k)

theta <- round(stats::rlnorm(n, 2))

alpha <- c(1, 2, 4, 1, 1)

custom_mmsbm <- mmsbm(
  theta = theta,
  B = B,
  alpha = alpha,
  expected_degree = 50
)

custom_mmsbm
#> Undirected Degree-Corrected Mixed Membership Stochastic Blockmodel
#> ------------------------------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 5
#> 
#> Traditional MMSBM parameterization:
#> 
#> Block memberships portions (Z): 100 x 5 [matrix] 
#> Degree heterogeneity (theta): 100 [numeric] 
#> Block propensities (alpha): 5 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 5 [dgeMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 5000
#> Expected degree: 50
#> Expected density: 1.0101

edgelist <- sample_edgelist(custom_mmsbm)
edgelist
#> # A tibble: 2,530 × 2
#>     from    to
#>    <int> <int>
#>  1    32    63
#>  2    11    19
#>  3    10    17
#>  4     1    15
#>  5     3    33
#>  6    13    14
#>  7     2    21
#>  8    14    46
#>  9     1     2
#> 10     2     8
#> # ℹ 2,520 more rows

# efficient eigendecompostion that leverages low-rank structure in
# E(A) so that you don't have to form E(A) to find eigenvectors,
# as E(A) is typically dense. computation is
# handled via RSpectra

population_eigs <- eigs_sym(custom_mmsbm)
svds(custom_mmsbm)$d
#> [1] 124.6072910   6.1422714   1.1009983   0.4416063   0.2927190
```
