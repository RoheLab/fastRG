# Create an undirected degree corrected stochastic blockmodel object

To specify a degree-corrected stochastic blockmodel, you must specify
the degree-heterogeneity parameters (via `n` or `theta`), the mixing
matrix (via `k` or `B`), and the relative block probabilities (optional,
via `pi`). We provide defaults for most of these options to enable rapid
exploration, or you can invest the effort for more control over the
model parameters. We **strongly recommend** setting the
`expected_degree` or `expected_density` argument to avoid large memory
allocations associated with sampling large, dense graphs.

## Usage

``` r
dcsbm(
  n = NULL,
  theta = NULL,
  k = NULL,
  B = NULL,
  ...,
  block_sizes = NULL,
  pi = NULL,
  sort_nodes = TRUE,
  force_identifiability = FALSE,
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
  will be symmetrized via the update `B := B + t(B)`. Defaults to
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

- block_sizes:

  (block sizes) Number of nodes in each block, as a vector of integers.
  Must match the dimensions of `B`, or `k` and must sum to `n`. Defaults
  to `NULL`, in which case blocks are made to be as balanced as
  possible. You can specify either `pi` or `block_sizes`, but not both.

- pi:

  (block sizes) Relative block probabilities. Must be positive, but do
  not need to sum to one, as they will be normalized internally. Must
  match the dimensions of `B` or `k`. Defaults to `NULL`, in which case
  the `block_sizes` argument will take precedence. Note that you can
  specify either `pi` or `block_sizes`, but should not specify both.

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block and by `theta`. Useful for plotting. Defaults to
  `TRUE`. When `TRUE`, nodes are first sorted by block membership, and
  then by degree-correction parameters within each block. Additionally,
  `pi` is sorted in increasing order, and the columns of the `B` matrix
  are permuted to match the new order of `pi`.

- force_identifiability:

  Logical indicating whether or not to normalize `theta` such that it
  sums to one within each block. Defaults to `FALSE`, since this
  behavior can be surprise when `theta` is set to a vector of all ones
  to recover the SBM case.

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

An `undirected_dcsbm` S3 object, a subclass of the
[`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)
with the following additional fields:

- `theta`: A numeric vector of degree-heterogeneity parameters.

- `z`: The community memberships of each node, as a
  [`factor()`](https://rdrr.io/r/base/factor.html). The factor will have
  `k` levels, where `k` is the number of communities in the stochastic
  blockmodel. There will not always necessarily be observed nodes in
  each community.

- `pi`: Sampling probabilities for each block.

- `sorted`: Logical indicating where nodes are arranged by block (and
  additionally by degree heterogeneity parameter) within each block.

## Generative Model

There are two levels of randomness in a degree-corrected stochastic
blockmodel. First, we randomly chose a block membership for each node in
the blockmodel. This is handled by `dcsbm()`. Then, given these block
memberships, we randomly sample edges between nodes. This second
operation is handled by
[`sample_edgelist()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md)
and
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md),
depending depending on your desired graph representation.

### Block memberships

Let \\z_i\\ represent the block membership of node \\i\\. To generate
\\z_i\\ we sample from a categorical distribution (note that this is a
special case of a multinomial) with parameter \\\pi\\, such that
\\\pi_i\\ represents the probability of ending up in the ith block.
Block memberships for each node are independent.

### Degree heterogeneity

In addition to block membership, the DCSBM also allows nodes to have
different propensities for edge formation. We represent this propensity
for node \\i\\ by a positive number \\\theta_i\\. Typically the
\\\theta_i\\ are constrained to sum to one for identifiability purposes,
but this doesn't really matter during sampling (i.e. without the sum
constraint scaling \\B\\ and \\\theta\\ has the same effect on edge
probabilities, but whether \\B\\ or \\\theta\\ is responsible for this
change is uncertain).

### Edge formulation

Once we know the block memberships \\z\\ and the degree heterogeneity
parameters \\theta\\, we need one more ingredient, which is the baseline
intensity of connections between nodes in block `i` and block `j`. Then
each edge \\A\_{i,j}\\ is Poisson distributed with parameter

\$\$ \lambda\[i, j\] = \theta_i \cdot B\_{z_i, z_j} \cdot \theta_j. \$\$

## See also

Other stochastic block models:
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

Other undirected graphs:
[`chung_lu()`](https://rohelab.github.io/fastRG/dev/reference/chung_lu.md),
[`erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/erdos_renyi.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

## Examples

``` r
set.seed(27)

lazy_dcsbm <- dcsbm(n = 100, k = 5, expected_density = 0.01)
#> Generating random degree heterogeneity parameters `theta` from a LogNormal(2, 1) distribution. This distribution may change in the future. Explicitly set `theta` for reproducible results.
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.
lazy_dcsbm
#> Undirected Degree-Corrected Stochastic Blockmodel
#> -------------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 5
#> 
#> Traditional DCSBM parameterization:
#> 
#> Block memberships (z): 100 [factor] 
#> Degree heterogeneity (theta): 100 [numeric] 
#> Block probabilities (pi): 5 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 5 [dgCMatrix] 
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

dense_lazy_dcsbm <- dcsbm(n = 50, k = 3, expected_density = 0.8)
#> Generating random degree heterogeneity parameters `theta` from a LogNormal(2, 1) distribution. This distribution may change in the future. Explicitly set `theta` for reproducible results.
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.
dense_lazy_dcsbm
#> Undirected Degree-Corrected Stochastic Blockmodel
#> -------------------------------------------------
#> 
#> Nodes (n): 50 (arranged by block)
#> Blocks (k): 3
#> 
#> Traditional DCSBM parameterization:
#> 
#> Block memberships (z): 50 [factor] 
#> Degree heterogeneity (theta): 50 [numeric] 
#> Block probabilities (pi): 3 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 50 x 3 [dgCMatrix] 
#> S: 3 x 3 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 980
#> Expected degree: 19.6
#> Expected density: 0.8

# explicitly setting the degree heterogeneity parameter,
# mixing matrix, and relative community sizes rather
# than using randomly generated defaults

k <- 5
n <- 100
B <- matrix(stats::runif(k * k), nrow = k, ncol = k)

theta <- round(stats::rlnorm(n, 2))

pi <- c(1, 2, 4, 1, 1)

custom_dcsbm <- dcsbm(
  theta = theta,
  B = B,
  pi = pi,
  expected_degree = 50
)

custom_dcsbm
#> Undirected Degree-Corrected Stochastic Blockmodel
#> -------------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 5
#> 
#> Traditional DCSBM parameterization:
#> 
#> Block memberships (z): 100 [factor] 
#> Degree heterogeneity (theta): 100 [numeric] 
#> Block probabilities (pi): 5 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 5 [dgCMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 5000
#> Expected degree: 50
#> Expected density: 1.0101

edgelist <- sample_edgelist(custom_dcsbm)
edgelist
#> # A tibble: 4,925 × 2
#>     from    to
#>    <int> <int>
#>  1     8    12
#>  2     3     7
#>  3     3     6
#>  4     3     6
#>  5     1    10
#>  6     2     9
#>  7     2     7
#>  8     1     4
#>  9     2     3
#> 10     1     5
#> # ℹ 4,915 more rows


dcsbm_explicit_block_sizes <- dcsbm(
  theta = rexp(100, 1 / 3) + 1,
  B = B,
  block_sizes = c(13, 17, 40, 14, 16),
  expected_degree = 5
)

# respects block sizes
summary(dcsbm_explicit_block_sizes$z)
#> block1 block2 block3 block4 block5 
#>     13     17     40     14     16 

# efficient eigendecompostion that leverages low-rank structure in
# E(A) so that you don't have to form E(A) to find eigenvectors,
# as E(A) is typically dense. computation is
# handled via RSpectra

population_eigs <- eigs_sym(custom_dcsbm)

```
