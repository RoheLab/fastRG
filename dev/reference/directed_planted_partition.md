# Create a directed planted partition object

To specify a directed planted partition model, you must specify the
number of nodes (via `n`), the number of blocks (via `k`), the edge
probabilities (either via `within_block` and `between_block`, or via `a`
and `b`), and the block sizes (optional, via `pi` or `block_sizes`). By
default, the communities a node sends edges from and receives edges into
are sampled independently; set `shared_blocks = TRUE` to give each node
a single community. We **strongly recommend** setting the
`expected_out_degree`, `expected_in_degree`, or `expected_density`
argument to avoid large memory allocations associated with sampling
large, dense graphs.

## Usage

``` r
directed_planted_partition(
  n,
  k,
  ...,
  within_block = NULL,
  between_block = NULL,
  a = NULL,
  b = NULL,
  block_sizes = NULL,
  pi = NULL,
  shared_blocks = FALSE,
  sort_nodes = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
)
```

## Arguments

- n:

  The number of nodes in the network. Must be a positive integer. This
  argument is required.

- k:

  Number of planted partitions, as a positive integer. There are `k`
  outgoing blocks and `k` incoming blocks. Must be less than or equal to
  `n`. This argument is required.

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

- within_block:

  Probability of a directed edge from a node in outgoing block `i` to a
  node in incoming block `i`, i.e. the diagonal of the mixing matrix
  `B`. Must be between zero and one. Must specify either `within_block`
  and `between_block`, or `a` and `b` to determine edge probabilities.

- between_block:

  Probability of a directed edge from a node in outgoing block `i` to a
  node in incoming block `j`, for `i != j`, i.e. the off-diagonal of the
  mixing matrix `B`. Must be between zero and one. Must specify either
  `within_block` and `between_block`, or `a` and `b` to determine edge
  probabilities.

- a:

  Number such that `a / n` is the probability of edges within a block.
  Useful for sparse graphs. Must specify either `within_block` and
  `between_block`, or `a` and `b` to determine edge probabilities.

- b:

  Number such that `b / n` is the probability of edges between blocks.
  Useful for sparse graphs. Must specify either `within_block` and
  `between_block`, or `a` and `b` to determine edge probabilities.

- block_sizes:

  (block sizes) Number of nodes in each block, as a vector of
  non-negative integers. Must have length `k` and must sum to `n`. Used
  for both the outgoing and the incoming blocks: exactly
  `block_sizes[j]` nodes are assigned to outgoing block `j`, and exactly
  `block_sizes[j]` nodes are assigned to incoming block `j`. Defaults to
  `NULL`, in which case blocks are made to be as balanced as possible
  (unless `pi` is specified). You can specify either `pi` or
  `block_sizes`, but not both.

- pi:

  (block sizes) Relative block probabilities, used for both the outgoing
  and the incoming blocks. Must be positive, but do not need to sum to
  one, as they will be normalized internally. Must have length `k`. When
  specified, block memberships are sampled independently for each node
  with these probabilities. Defaults to `NULL`, in which case the
  `block_sizes` argument (or its balanced default) takes precedence. You
  can specify either `pi` or `block_sizes`, but not both.

- shared_blocks:

  Logical indicating whether each node should belong to the same
  outgoing and incoming block. Defaults to `FALSE`, in which case
  outgoing and incoming block memberships are sampled independently.
  When `TRUE`, a single block membership vector is sampled and used for
  both, so that `z_out` and `z_in` are identical. See
  [`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md)
  for details.

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block. Useful for plotting. Defaults to `TRUE`. When
  `TRUE`, blocks are relabeled in increasing order of `pi` (or
  `block_sizes`), and the rows of the expected adjacency matrix are
  grouped by outgoing block while the columns are grouped by incoming
  block. When `shared_blocks = TRUE`, rows and columns are permuted
  together, so that row `i` and column `i` refer to the same node. When
  `shared_blocks = FALSE` (the default), rows and columns are sorted
  *independently*, so row `i` and column `i` are generally not the same
  node: the diagonal of `A` (self-loops) is then not meaningful as "node
  `i` to itself", and `z_out[i]` and `z_in[i]` describe different nodes.
  Use `shared_blocks = TRUE` or `sort_nodes = FALSE` if you need rows
  and columns to index the same nodes. See
  [`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md)
  for details.

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

A `directed_planted_partition` S3 object, which is a subclass of the
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md)
object (with `theta_out` and `theta_in` identically one), with
additional fields:

- `within_block`: The probability of edge formation within a block, as
  specified (or as computed from `a / n`), before any rescaling by
  `expected_out_degree`, `expected_in_degree` or `expected_density`.

- `between_block`: The probability of edge formation between two
  distinct blocks, as specified (or as computed from `b / n`), before
  any rescaling by `expected_out_degree`, `expected_in_degree` or
  `expected_density`.

## Details

A directed planted partition model is a directed stochastic blockmodel
with `k` outgoing blocks and `k` incoming blocks, in which the diagonal
and the off-diagonal of the mixing matrix `B` are both constant. This
means that the probability of an edge from node \\i\\ to node \\j\\
depends only on whether the outgoing block of \\i\\ matches the incoming
block of \\j\\, but the particular blocks themselves don't have any
impact apart from this. It is a special case of
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md)
where all degree-heterogeneity parameters are equal to one.

## Generative Model

### Block memberships

Each node \\i\\ has an outgoing block membership \\z\_{out, i}\\, which
governs the edges it sends, and an incoming block membership \\z\_{in,
i}\\, which governs the edges it receives. Both take values in
`1, ..., k`. By default, the outgoing memberships are a uniformly random
assignment of nodes to blocks with block sizes that are as balanced as
possible (as in
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md)),
or with exactly `block_sizes[j]` nodes in block `j` when `block_sizes`
is specified. The incoming memberships are assigned in the same way,
independently of the outgoing memberships. When `pi` is specified,
memberships are instead sampled independently across nodes from a
categorical distribution with parameter \\\pi\\ (again, independently
for outgoing and incoming blocks).

When `shared_blocks = TRUE`, a single membership vector is sampled and
\\z\_{out, i} = z\_{in, i}\\ for every node.

### Mixing matrix

The mixing matrix is the `k` by `k` matrix with `within_block` (or
`a / n`) on the diagonal and `between_block` (or `b / n`) everywhere
else, i.e.

\$\$ B = p\_{between} \mathbf{1} \mathbf{1}^T + (p\_{within} -
p\_{between}) I_k, \$\$

where \\p\_{within}\\ is `within_block` and \\p\_{between}\\ is
`between_block`. If you specify `expected_out_degree`,
`expected_in_degree` or `expected_density`, `B` is rescaled by a
constant, so only the ratio of `within_block` to `between_block`
matters.

### Edge formation

Given the block memberships, each directed edge \\A\_{ij}\\ from node
\\i\\ to node \\j\\ is sampled independently (for every ordered pair
\\(i, j)\\, so the adjacency matrix is not symmetric) with expectation

\$\$ E\[A\_{ij}\] = B\_{z\_{out, i}, z\_{in, j}}. \$\$

Edges are Poisson distributed when `poisson_edges = TRUE` (yielding a
multigraph) and Bernoulli distributed with the same probabilities when
`poisson_edges = FALSE`. Self-loops are allowed or not according to
`allow_self_loops`. Note that with the default `shared_blocks = FALSE`
and `sort_nodes = TRUE`, rows and columns of \\A\\ are sorted
independently, so the diagonal entry \\A\_{ii}\\ pairs the \\i\\-th
sender with the \\i\\-th receiver rather than node \\i\\ with itself
(see `sort_nodes`).

### Interpretation

When `within_block > between_block`, nodes in outgoing block \\s\\
preferentially send edges to nodes in incoming block \\s\\. With the
default independent memberships, this is a co-clustering structure: the
communities that describe who a node sends edges to are unrelated to the
communities that describe who sends edges to it. With
`shared_blocks = TRUE`, each node has one community, and edges in both
directions are more common within communities than between them.

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

Other directed graphs:
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md),
[`directed_erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/directed_erdos_renyi.md)

## Examples

``` r

set.seed(27)

lazy_dpp <- directed_planted_partition(
  n = 100,
  k = 5,
  expected_density = 0.01,
  within_block = 0.1,
  between_block = 0.01
)

lazy_dpp
#> Directed Planted Partition Model
#> --------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 5
#> Shared blocks: FALSE (z_out and z_in are independent)
#> 
#> Traditional planted partition parameterization:
#> 
#> Within block edge probability: 0.1 
#> Between block edge probability: 0.01 
#> Block memberships (z_out): 100 [factor] 
#> Block memberships (z_in): 100 [factor] 
#> Block probabilities (pi_out): 5 [numeric] 
#> Block probabilities (pi_in): 5 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 5 [dgCMatrix] 
#> S: 5 x 5 [dsyMatrix] 
#> Y: 100 x 5 [dgCMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 100
#> Expected out degree: 1
#> Expected in degree: 1
#> Expected density: 0.01

# sparse parameterization, with Bernoulli edges

sparse_dpp <- directed_planted_partition(
  n = 1000,
  k = 3,
  a = 10,
  b = 4,
  poisson_edges = FALSE
)

A <- sample_sparse(sparse_dpp)

# only zeroes and ones!
range(A)
#> [1] 0 1

# each node belongs to a single community, which determines both
# the edges it sends and the edges it receives

shared_dpp <- directed_planted_partition(
  n = 100,
  k = 4,
  within_block = 0.3,
  between_block = 0.02,
  shared_blocks = TRUE
)

identical(shared_dpp$z_out, shared_dpp$z_in)
#> [1] TRUE

# exact block sizes

sized_dpp <- directed_planted_partition(
  n = 90,
  k = 3,
  within_block = 0.3,
  between_block = 0.05,
  block_sizes = c(20, 30, 40)
)

table(sized_dpp$z_out)
#> 
#>  1  2  3 
#> 20 30 40 
table(sized_dpp$z_in)
#> 
#>  1  2  3 
#> 20 30 40 

population_svd <- svds(sized_dpp)
population_svd$d
#> [1] 13.239869  8.300499  5.459632
```
