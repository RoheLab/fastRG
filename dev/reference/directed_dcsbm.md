# Create a directed degree corrected stochastic blockmodel object

To specify a degree-corrected stochastic blockmodel, you must specify
the degree-heterogeneity parameters (via `n` or `theta_out` and
`theta_in`), the mixing matrix (via `k_out` and `k_in`, or `B`), and the
block sizes (optional, either as relative block probabilities via
`pi_out` and `pi_in`, or as exact block counts via `block_sizes_out` and
`block_sizes_in`). By default, outgoing and incoming block memberships
and degree parameters are unrelated (a co-factor model); use
`shared_blocks` and `shared_theta` to tie them together. We provide
defaults for most of these options to enable rapid exploration, or you
can invest the effort for more control over the model parameters. We
**strongly recommend** setting the `expected_out_degree`,
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
  block_sizes_out = NULL,
  block_sizes_in = NULL,
  pi_out = NULL,
  pi_in = NULL,
  shared_blocks = FALSE,
  shared_theta = FALSE,
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
  `theta_in` together, but not both. When `shared_theta = TRUE`, a
  single vector is generated and used for both `theta_out` and
  `theta_in`.

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
  once, in which case the number of outgoing blocks is taken from `B`.

- k_in:

  (mixing matrix) The number of incoming blocks in the blockmodel. Use
  when you don't want to specify the mixing-matrix by hand. When `k_in`
  is specified, the elements of `B` are drawn randomly from a
  `Uniform(0, 1)` distribution. This is subject to change, and may not
  be reproducible. `k_in` defaults to `NULL`. You may specify all three
  at once, in which case the number of incoming blocks is taken from
  `B`.

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

- block_sizes_out:

  (block sizes) Number of nodes in each outgoing block, as a vector of
  non-negative integers. Must match the number of rows of `B` (or
  `k_out`) and must sum to `n` (or `length(theta_out)`). When specified,
  exactly `block_sizes_out[j]` nodes are assigned to outgoing block `j`,
  and which nodes land in which block is uniformly random. Defaults to
  `NULL`, in which case outgoing blocks are made to be as balanced as
  possible, unless `pi_out` is specified. You can specify either
  `pi_out` or `block_sizes_out`, but not both.

- block_sizes_in:

  (block sizes) Number of nodes in each incoming block, as a vector of
  non-negative integers. Must match the number of columns of `B` (or
  `k_in`) and must sum to `n` (or `length(theta_in)`). When specified,
  exactly `block_sizes_in[j]` nodes are assigned to incoming block `j`,
  and which nodes land in which block is uniformly random. Defaults to
  `NULL`, in which case incoming blocks are made to be as balanced as
  possible, unless `pi_in` is specified. You can specify either `pi_in`
  or `block_sizes_in`, but not both.

- pi_out:

  (block sizes) Relative block probabilities for outgoing blocks. Must
  be positive, but do not need to sum to one, as they will be normalized
  internally. Must match the rows of `B`, or `k_out`. When specified,
  outgoing block memberships are sampled independently for each node
  with these probabilities. Defaults to `NULL`, in which case the
  `block_sizes_out` argument (or its balanced default) takes precedence.
  You can specify either `pi_out` or `block_sizes_out`, but not both.

- pi_in:

  (block sizes) Relative block probabilities for incoming blocks. Must
  be positive, but do not need to sum to one, as they will be normalized
  internally. Must match the columns of `B`, or `k_in`. When specified,
  incoming block memberships are sampled independently for each node
  with these probabilities. Defaults to `NULL`, in which case the
  `block_sizes_in` argument (or its balanced default) takes precedence.
  You can specify either `pi_in` or `block_sizes_in`, but not both.

- shared_blocks:

  (block memberships) Logical indicating whether each node should belong
  to the same outgoing and incoming block. Defaults to `FALSE`, in which
  case outgoing and incoming block memberships are sampled
  independently. When `TRUE`, a single block membership vector is
  sampled and used for both, so that `z_out` and `z_in` are identical.
  This requires the same number of outgoing and incoming blocks, i.e.
  `B` must be square (or `k_out` must equal `k_in`; specifying only one
  of them is enough). It is enough to specify only one of `pi_out` and
  `pi_in` (or of `block_sizes_out` and `block_sizes_in`); if both are
  specified, they must agree.

- shared_theta:

  (degree heterogeneity) Logical indicating whether each node should
  have the same propensity to form outgoing and incoming edges. Defaults
  to `FALSE`, in which case `theta_out` and `theta_in` are unrelated.
  When `TRUE`, a single degree-heterogeneity vector is used for both
  `theta_out` and `theta_in`: either specify `n` (one vector is then
  generated at random), or only one of `theta_out` and `theta_in`. If
  you specify both `theta_out` and `theta_in`, they must be identical.
  Requires `shared_blocks = TRUE`, so that the stored `theta_out` and
  `theta_in` are identical (also after normalization by
  `force_identifiability`).

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block. Useful for plotting. Defaults to `TRUE`. When
  `TRUE`, rows of the expected adjacency matrix are first sorted by
  outgoing block membership, and then by outgoing degree-correction
  parameters within each outgoing block. A similar sorting procedure
  occurs independently for the columns, according to the incoming
  blocks. Additionally, `pi_out` and `pi_in` are sorted in increasing
  order, and the rows and columns of the `B` matrix are permuted to
  match the new orderings. When `block_sizes_out` or `block_sizes_in`
  are specified, they are permuted in the same way, so each block keeps
  exactly the requested number of nodes and its row (or column) of `B`.

  **Rows and columns index different nodes.** In the default co-factor
  model (`shared_blocks = FALSE`), the rows and the columns are permuted
  *independently*, so after sorting, row `i` and column `i` of the
  adjacency matrix are generally not the same node. In particular, the
  diagonal of `A` (self-loops, see `allow_self_loops`) does not mean
  "node `i` to itself", and `X[i, ]` and `Y[i, ]` describe different
  nodes. If you need rows and columns to index the same nodes, use
  `shared_blocks = TRUE` or `sort_nodes = FALSE`.

  When `shared_blocks = TRUE`, rows and columns are instead permuted
  together, so that row `i` and column `i` of the expected adjacency
  matrix refer to the same node: nodes are sorted by block, and then by
  decreasing `theta_out` and `theta_in` within each block.

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

- `theta_out`: A numeric vector of outgoing degree-heterogeneity
  parameters.

- `theta_in`: A numeric vector of incoming degree-heterogeneity
  parameters.

- `z_out`: The outgoing community memberships of each node, as a
  [`factor()`](https://rdrr.io/r/base/factor.html). The factor will have
  `k_out` levels, where `k_out` is the number of outgoing communities in
  the stochastic blockmodel. There will not always necessarily be
  observed nodes in each community.

- `z_in`: The incoming community memberships of each node, as a
  [`factor()`](https://rdrr.io/r/base/factor.html). The factor will have
  `k_in` levels, where `k_in` is the number of incoming communities in
  the stochastic blockmodel. There will not always necessarily be
  observed nodes in each community.

- `pi_out`: Sampling probabilities for each outgoing community.

- `pi_in`: Sampling probabilities for each incoming community.

- `sorted`: Logical indicating where nodes are arranged by block (and
  additionally by degree heterogeneity parameter) within each block.

- `shared_blocks`: Logical indicating whether outgoing and incoming
  block memberships are shared, in which case `z_out` and `z_in` are
  identical.

- `shared_theta`: Logical indicating whether outgoing and incoming
  degree-heterogeneity parameters are shared, in which case `theta_out`
  and `theta_in` are identical.

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
sample from a categorical distribution with parameter \\\pi\_{out}\\. To
generate \\y\\ we sample from a categorical distribution with parameter
\\\pi\_{in}\\. Block memberships are independent across nodes. Incoming
and outgoing block memberships of the same node are also independent.

These categorical draws are used when `pi_out` (`pi_in`) is specified.
Otherwise, the outgoing (incoming) block memberships are a uniformly
random assignment of nodes to blocks with fixed block sizes: either
exactly the sizes given in `block_sizes_out` (`block_sizes_in`), or, by
default, block sizes that are as balanced as possible (every block gets
`floor(n / k_out)` nodes, and the remainder is spread over the first
blocks, one node each), as in
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md).
Fixed block sizes are useful, for instance, to construct blockmodels
with exactly repeated singular values.

### Degree heterogeneity

In addition to block membership, the DCSBM also nodes to have different
propensities for incoming and outgoing edge formation. We represent the
propensity to form incoming edges for a given node by a positive number
\\\theta\_{out}\\. We represent the propensity to form outgoing edges
for a given node by a positive number \\\theta\_{in}\\. Typically the
\\\theta\_{out}\\ (and \\\theta\_{in}\\) across all nodes are
constrained to sum to one for identifiability purposes, but this doesn't
really matter during sampling.

### Edge formulation

Once we know the block memberships \\x\\ and \\y\\ and the degree
heterogeneity parameters \\\theta\_{in}\\ and \\\theta\_{out}\\, we need
one more ingredient, which is the baseline intensity of connections
between nodes in block `i` and block `j`. Then each edge forms
independently according to a Poisson distribution with parameters

\$\$ \lambda = \theta\_{in} \cdot B\_{x, y} \cdot \theta\_{out}. \$\$

## Co-factor models versus shared blocks and degree heterogeneity

By default (`shared_blocks = FALSE` and `shared_theta = FALSE`),
`directed_dcsbm()` is a *co-factor* model: every node has an outgoing
block \\z\_{out}\\ and degree parameter \\\theta\_{out}\\ that govern
the edges it sends (the rows of \\A\\), and a separate incoming block
\\z\_{in}\\ and degree parameter \\\theta\_{in}\\ that govern the edges
it receives (the columns of \\A\\). These are drawn independently, so
the communities a node sends edges from and receives edges into are
unrelated, and the two sets of communities can even have different
numbers of blocks. This is the model targeted by co-clustering methods
that estimate separate sending and receiving communities.

In the co-factor model, `sort_nodes = TRUE` (the default) sorts the rows
and the columns of \\A\\ *independently*, so row \\i\\ and column \\i\\
are generally not the same node. Consequently the diagonal of \\A\\, and
so the self-loops controlled by `allow_self_loops`, is not meaningful as
"node \\i\\ to itself": entry \\A\_{ii}\\ connects the \\i\\-th sender
to the \\i\\-th receiver in their respective orderings. If you need rows
and columns to index the same nodes (for example, to interpret
self-loops or to read the graph as a directed graph on a single set of
nodes), use `shared_blocks = TRUE` or `sort_nodes = FALSE`.

Often it is more natural to think of each node as belonging to a single
community, which shapes both the edges it sends and the edges it
receives. Setting `shared_blocks = TRUE` draws one membership vector
\\z\\ with \\z\_{out} = z\_{in} = z\\, so that \\E\[A\_{ij}\] =
\theta\_{out, i} B\_{z_i, z_j} \theta\_{in, j}\\, where an asymmetric
`B` encodes directed relationships between the communities. With shared
blocks you can additionally set `shared_theta = TRUE`, which uses one
degree-heterogeneity vector \\\theta\\ for both \\\theta\_{out}\\ and
\\\theta\_{in}\\, so nodes that send many edges also tend to receive
many edges. (`shared_theta = TRUE` requires `shared_blocks = TRUE`.)
With shared blocks, `sort_nodes = TRUE` permutes rows and columns
together, so that row `i` and column `i` of \\A\\ always refer to the
same node.

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`directed_planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/directed_planted_partition.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

Other directed graphs:
[`directed_erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/directed_erdos_renyi.md),
[`directed_planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/directed_planted_partition.md)

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
#> Outgoing Blocks (k_out): 5
#> Incoming Blocks (k_in): 8
#> Shared blocks: FALSE (z_out and z_in are independent)
#> Shared degree heterogeneity: FALSE (theta_out and theta_in are separate)
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
#> Expected out degree: 1
#> Expected in degree: 1
#> Expected density: 0.01

population_svd <- svds(ddcsbm)

# specify the exact number of nodes in each outgoing and incoming block,
# rather than sampling block memberships

ddcsbm_explicit_block_sizes <- directed_dcsbm(
  theta_out = rexp(100, 1 / 3) + 1,
  theta_in = rexp(100, 1 / 3) + 1,
  B = B,
  block_sizes_out = c(10, 15, 20, 25, 30),
  block_sizes_in = c(5, 5, 10, 10, 15, 15, 20, 20),
  expected_density = 0.01
)

# respects block sizes
table(ddcsbm_explicit_block_sizes$z_out)
#> 
#>  1  2  3  4  5 
#> 10 15 20 25 30 
table(ddcsbm_explicit_block_sizes$z_in)
#> 
#>  1  2  3  4  5  6  7  8 
#>  5  5 10 10 15 15 20 20 

# each node belongs to the same outgoing and incoming community, and has
# the same propensity to send and receive edges

B_square <- matrix(0.1, nrow = 4, ncol = 4)
diag(B_square) <- 0.8
B_square[1, 2] <- 0.5 # block 1 sends many edges to block 2

shared_ddcsbm <- directed_dcsbm(
  theta_out = rexp(100, 1 / 3) + 1,
  B = B_square,
  shared_blocks = TRUE,
  shared_theta = TRUE,
  expected_density = 0.05
)

shared_ddcsbm
#> Directed Degree-Corrected Stochastic Blockmodel
#> -----------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Outgoing Blocks (k_out): 4
#> Incoming Blocks (k_in): 4
#> Shared blocks: TRUE (z_out and z_in are identical)
#> Shared degree heterogeneity: TRUE (theta_out and theta_in are identical)
#> 
#> Traditional DCSBM parameterization:
#> 
#> Block memberships (z_out): 100 [factor] 
#> Block memberships (z_in): 100 [factor] 
#> Degree heterogeneity (theta_out): 100 [numeric] 
#> Degree heterogeneity (theta_in): 100 [numeric] 
#> Block probabilities (pi_out): 4 [numeric] 
#> Block probabilities (pi_in): 4 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 4 [dgCMatrix] 
#> S: 4 x 4 [dgeMatrix] 
#> Y: 100 x 4 [dgCMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 500
#> Expected out degree: 5
#> Expected in degree: 5
#> Expected density: 0.05

identical(shared_ddcsbm$z_out, shared_ddcsbm$z_in)
#> [1] TRUE
```
