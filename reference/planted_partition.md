# Create an undirected planted partition object

To specify a planted partition model, you must specify the number of
nodes (via `n`), the mixing matrix (optional, either via
`within_block/between_block` or `a/b`), and the relative block
probabilites (optional, via `pi`). We provide defaults for most of these
options to enable rapid exploration, or you can invest the effort for
more control over the model parameters. We **strongly recommend**
setting the `expected_degree` or `expected_density` argument to avoid
large memory allocations associated with sampling large, dense graphs.

## Usage

``` r
planted_partition(
  n,
  k,
  ...,
  within_block = NULL,
  between_block = NULL,
  a = NULL,
  b = NULL,
  block_sizes = NULL,
  pi = NULL,
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

  Number of planted partitions, as a positive integer. This argument is
  required.

- ...:

  Arguments passed on to
  [`undirected_factor_model`](https://rohelab.github.io/fastRG/reference/undirected_factor_model.md)

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

- within_block:

  Probability of within block edges. Must be strictly between zero and
  one. Must specify either `within_block` and `between_block`, or `a`
  and `b` to determine edge probabilities.

- between_block:

  Probability of between block edges. Must be strictly between zero and
  one. Must specify either `within_block` and `between_block`, or `a`
  and `b` to determine edge probabilities.

- a:

  Integer such that `a/n` is the probability of edges within a block.
  Useful for sparse graphs. Must specify either `within_block` and
  `between_block`, or `a` and `b` to determine edge probabilities.

- b:

  Integer such that `b/n` is the probability of edges between blocks.
  Useful for sparse graphs. Must specify either `within_block` and
  `between_block`, or `a` and `b` to determine edge probabilities.

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

An `undirected_planted_partition` S3 object, which is a subclass of the
[`sbm()`](https://rohelab.github.io/fastRG/reference/sbm.md) object,
with additional fields:

- `within_block`: The probability of edge formation within a block.

- `between_block`: The probability of edge formation between two
  distinct blocks.

## Details

A planted partition model is stochastic blockmodel in which the diagonal
and the off-diagonal of the mixing matrix `B` are both constant. This
means that edge probabilities depend only on whether two nodes belong to
the same block, or to different blocks, but the particular blocks
themselves don't have any impact apart from this.

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/reference/dcsbm.md),
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/reference/directed_dcsbm.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/reference/overlapping_sbm.md),
[`sbm()`](https://rohelab.github.io/fastRG/reference/sbm.md)

Other undirected graphs:
[`chung_lu()`](https://rohelab.github.io/fastRG/reference/chung_lu.md),
[`dcsbm()`](https://rohelab.github.io/fastRG/reference/dcsbm.md),
[`erdos_renyi()`](https://rohelab.github.io/fastRG/reference/erdos_renyi.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/reference/overlapping_sbm.md),
[`sbm()`](https://rohelab.github.io/fastRG/reference/sbm.md)

## Examples

``` r
set.seed(27)

lazy_pp <- planted_partition(
  n = 100,
  k = 5,
  expected_density = 0.01,
  within_block = 0.1,
  between_block = 0.01
)

lazy_pp
#> Undirected Stochastic Blockmodel
#> --------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 5
#> 
#> Traditional SBM parameterization:
#> 
#> Block memberships (z): 100 [factor] 
#> Block probabilities (pi): 5 [numeric] 
#> Factor model parameterization:
#> 
#> X: 100 x 5 [dgCMatrix] 
#> S: 5 x 5 [dsyMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 50
#> Expected degree: 0.5
#> Expected density: 0.01
```
