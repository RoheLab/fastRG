# Create a directed degree corrected stochastic blockmodel object

To specify a degree-corrected stochastic blockmodel, you must specify
the degree-heterogeneity parameters (via `n` or `theta_out` and
`theta_in`), the mixing matrix (via `k_out` and `k_in`, or `B`), and the
relative block probabilities (optional, via `p_out` and `pi_in`). We
provide defaults for most of these options to enable rapid exploration,
or you can invest the effort for more control over the model parameters.
We **strongly recommend** setting the `expected_out_degree`,
`expected_in_degree`, or `expected_density` argument to avoid large
memory allocations associated with sampling large, dense graphs.

## Usage

``` r
directed_dcsbm(
  n = NULL,
  theta_out = NULL,
  theta_in = NULL,
  k_out = NULL,
  k_in = NULL,
  B = NULL,
  ...,
  pi_out = rep(1/k_out, k_out),
  pi_in = rep(1/k_in, k_in),
  sort_nodes = TRUE,
  force_identifiability = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
)
```

## Arguments

- n:

  (degree heterogeneity) The number of nodes in the blockmodel. Use when
  you don't want to specify the degree-heterogeneity parameters
  `theta_out` and `theta_in` by hand. When `n` is specified, `theta_out`
  and `theta_in` are randomly generated from a `LogNormal(2, 1)`
  distribution. This is subject to change, and may not be reproducible.
  `n` defaults to `NULL`. You must specify either `n` or `theta_out` and
  `theta_in` together, but not both.

- theta_out:

  (degree heterogeneity) A numeric vector explicitly specifying the
  degree heterogeneity parameters. This implicitly determines the number
  of nodes in the resulting graph, i.e. it will have `length(theta_out)`
  nodes. Must be positive. Setting to a vector of ones recovers a
  stochastic blockmodel without degree correction. Defaults to `NULL`.
  You must specify either `n` or `theta_out` and `theta_in` together,
  but not both. `theta_out` controls outgoing degree propensity, or,
  equivalently, row sums of the adjacency matrix.

- theta_in:

  (degree heterogeneity) A numeric vector explicitly specifying the
  degree heterogeneity parameters. This implicitly determines the number
  of nodes in the resulting graph, i.e. it will have `length(theta)`
  nodes. Must be positive. Setting to a vector of ones recovers a
  stochastic blockmodel without degree correction. Defaults to `NULL`.
  You must specify either `n` or `theta_out` and `theta_in` together,
  but not both. `theta_in` controls incoming degree propensity, or,
  equivalently, column sums of the adjacency matrix.

- k_out:

  (mixing matrix) The number of outgoing blocks in the blockmodel. Use
  when you don't want to specify the mixing-matrix by hand. When `k_out`
  is specified, the elements of `B` are drawn randomly from a
  `Uniform(0, 1)` distribution. This is subject to change, and may not
  be reproducible. `k_out` defaults to `NULL`. You must specify either
  `k_out` and `k_in` together, or `B`. You may specify all three at
  once, in which case `k_out` is only used to set `pi_out` (when
  `pi_out` is left at its default argument value).

- k_in:

  (mixing matrix) The number of incoming blocks in the blockmodel. Use
  when you don't want to specify the mixing-matrix by hand. When `k_in`
  is specified, the elements of `B` are drawn randomly from a
  `Uniform(0, 1)` distribution. This is subject to change, and may not
  be reproducible. `k_in` defaults to `NULL`. You may specify all three
  at once, in which case `k_in` is only used to set `pi_in` (when
  `pi_in` is left at its default argument value).

- B:

  (mixing matrix) A `k_out` by `k_in` matrix of block connection
  probabilities. The probability that a node in block `i` connects to a
  node in community `j` is `Poisson(B[i, j])`. `matrix` and `Matrix`
  objects are both acceptable. Defaults to `NULL`. You must specify
  either `k_out` and `k_in` together, or `B`, but not both.

- ...:

  Arguments passed on to
  [`directed_factor_model`](https://rohelab.github.io/fastRG/dev/reference/directed_factor_model.md)

  `expected_in_degree`

  :   If specified, the desired expected in degree of the graph.
      Specifying `expected_in_degree` simply rescales `S` to achieve
      this. Defaults to `NULL`. Specify only one of
      `expected_in_degree`, `expected_out_degree`, and
      `expected_density`.

  `expected_out_degree`

  :   If specified, the desired expected out degree of the graph.
      Specifying `expected_out_degree` simply rescales `S` to achieve
      this. Defaults to `NULL`. Specify only one of
      `expected_in_degree`, `expected_out_degree`, and
      `expected_density`.

  `expected_density`

  :   If specified, the desired expected density of the graph.
      Specifying `expected_density` simply rescales `S` to achieve this.
      Defaults to `NULL`. Specify only one of `expected_in_degree`,
      `expected_out_degree`, and `expected_density`.

- pi_out:

  (relative block probabilities) Relative block probabilities. Must be
  positive, but do not need to sum to one, as they will be normalized
  internally. Must match the rows of `B`, or `k_out`. Defaults to
  `rep(1 / k_out, k_out)`, or a balanced outgoing blocks.

- pi_in:

  (relative block probabilities) Relative block probabilities. Must be
  positive, but do not need to sum to one, as they will be normalized
  internally. Must match the columns of `B`, or `k_in`. Defaults to
  `rep(1 / k_in, k_in)`, or a balanced incoming blocks.

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block. Useful for plotting. Defaults to `TRUE`. When
  `TRUE`, rows of the expected adjacency matrix are first sorted by
  outgoing block membership, and then by incoming degree-correction
  parameters within each incoming block. A similar sorting procedure
  occurs independently from the columns, according to the incoming
  blocks. Additionally, `pi_out` and `pi_in` are sorted in increasing
  order, and the columns of the `B` matrix are permuted to match the new
  orderings.

- force_identifiability:

  Logical indicating whether or not to normalize `theta_out` such that
  it sums to one within each incoming block and `theta_in` such that it
  sums to one within each outgoing block. Defaults to `TRUE`.

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

A `directed_dcsbm` S3 object, a subclass of the
[`directed_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/directed_factor_model.md)
with the following additional fields:

- `theta_out`: A numeric vector of incoming community
  degree-heterogeneity parameters.

- `theta_in`: A numeric vector of outgoing community
  degree-heterogeneity parameters.

- `z_out`: The incoming community memberships of each node, as a
  [`factor()`](https://rdrr.io/r/base/factor.html). The factor will have
  `k_out` levels, where `k_out` is the number of incoming communities in
  the stochastic blockmodel. There will not always necessarily be
  observed nodes in each community.

- `z_in`: The outgoing community memberships of each node, as a
  [`factor()`](https://rdrr.io/r/base/factor.html). The factor will have
  `k_in` levels, where `k_in` is the number of outgoing communities in
  the stochastic blockmodel. There will not always necessarily be
  observed nodes in each community.

- `pi_out`: Sampling probabilities for each incoming community.

- `pi_in`: Sampling probabilities for each outgoing community.

- `sorted`: Logical indicating where nodes are arranged by block (and
  additionally by degree heterogeneity parameter) within each block.

## Generative Model

There are two levels of randomness in a directed degree-corrected
stochastic blockmodel. First, we randomly chose a incoming block
membership and an outgoing block membership for each node in the
blockmodel. This is handled by `directed_dcsbm()`. Then, given these
block memberships, we randomly sample edges between nodes. This second
operation is handled by
[`sample_edgelist()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md)
and
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md),
depending on your desired graph representation.

### Block memberships

Let \\x\\ represent the incoming block membership of a node and \\y\\
represent the outgoing block membership of a node. To generate \\x\\ we
sample from a categorical distribution with parameter \\\pi_out\\. To
generate \\y\\ we sample from a categorical distribution with parameter
\\\pi_in\\. Block memberships are independent across nodes. Incoming and
outgoing block memberships of the same node are also independent.

### Degree heterogeneity

In addition to block membership, the DCSBM also nodes to have different
propensities for incoming and outgoing edge formation. We represent the
propensity to form incoming edges for a given node by a positive number
\\\theta_out\\. We represent the propensity to form outgoing edges for a
given node by a positive number \\\theta_in\\. Typically the
\\\theta_out\\ (and \\theta_in\\) across all nodes are constrained to
sum to one for identifiability purposes, but this doesn't really matter
during sampling.

### Edge formulation

Once we know the block memberships \\x\\ and \\y\\ and the degree
heterogeneity parameters \\\theta\_{in}\\ and \\\theta\_{out}\\, we need
one more ingredient, which is the baseline intensity of connections
between nodes in block `i` and block `j`. Then each edge forms
independently according to a Poisson distribution with parameters

\$\$ \lambda = \theta\_{in} \* B\_{x, y} \* \theta\_{out}. \$\$

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

Other directed graphs:
[`directed_erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/directed_erdos_renyi.md)

## Examples

``` r
set.seed(27)

B <- matrix(0.2, nrow = 5, ncol = 8)
diag(B) <- 0.9

ddcsbm <- directed_dcsbm(
  n = 100,
  B = B,
  k_out = 5,
  k_in = 8,
  expected_density = 0.01
)
#> Generating random degree heterogeneity parameters `theta_out` and `theta_in` from LogNormal(2, 1) distributions. This distribution may change in the future. Explicitly set `theta_out` and `theta_in` for reproducible results.

ddcsbm
#> Directed Degree-Corrected Stochastic Blockmodel
#> -----------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Incoming Blocks (k_out): 5
#> Outgoing Blocks (k_in): 8
#> 
#> Traditional DCSBM parameterization:
#> 
#> Block memberships (z_out): 100 [factor] 
#> Block memberships (z_in): 100 [factor] 
#> Degree heterogeneity (theta_out): 100 [numeric] 
#> Degree heterogeneity (theta_in): 100 [numeric] 
#> Block probabilities (pi_out): 5 [numeric] 
#> Block probabilities (pi_in): 8 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 5 [dgCMatrix] 
#> S: 5 x 8 [dgeMatrix] 
#> Y: 100 x 8 [dgCMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 100
#> Expected in degree: 1
#> Expected out degree: 1
#> Expected density: 0.01

population_svd <- svds(ddcsbm)
```
