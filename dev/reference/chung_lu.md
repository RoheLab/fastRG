# Create an undirected Chung-Lu object

To specify a Chung-Lu graph, you must specify the degree-heterogeneity
parameters (via `n` or `theta`). We provide reasonable defaults to
enable rapid exploration or you can invest the effort for more control
over the model parameters. We **strongly recommend** setting the
`expected_degree` or `expected_density` argument to avoid large memory
allocations associated with sampling large, dense graphs.

## Usage

``` r
chung_lu(
  n = NULL,
  theta = NULL,
  ...,
  sort_nodes = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE,
  force_identifiability = FALSE
)
```

## Arguments

- n:

  (degree heterogeneity) The number of nodes in the graph. Use when you
  don't want to specify the degree-heterogeneity parameters `theta` by
  hand. When `n` is specified, `theta` is randomly generated from a
  `LogNormal(2, 1)` distribution. This is subject to change, and may not
  be reproducible. `n` defaults to `NULL`. You must specify either `n`
  or `theta`, but not both.

- theta:

  (degree heterogeneity) A numeric vector explicitly specifying the
  degree heterogeneity parameters. This implicitly determines the number
  of nodes in the resulting graph, i.e. it will have `length(theta)`
  nodes. Must be positive. Setting to a vector of ones recovers an erdos
  renyi graph. Defaults to `NULL`. You must specify either `n` or
  `theta`, but not both.

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

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block and by `theta`. Useful for plotting. Defaults to
  `TRUE`.

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

- force_identifiability:

  Logical indicating whether or not to normalize `theta` such that it
  sums to one within each block. Defaults to `FALSE`, since this
  behavior can be surprise when `theta` is set to a vector of all ones
  to recover the SBM case.

## Value

An `undirected_chung_lu` S3 object, a subclass of
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md).

## See also

Other undirected graphs:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/erdos_renyi.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)

## Examples

``` r
set.seed(27)

cl <- chung_lu(n = 100, expected_density = 0.01)
#> Generating random degree heterogeneity parameters `theta` from a LogNormal(2, 1) distribution. This distribution may change in the future. Explicitly set `theta` for reproducible results.
cl
#> Undirected Degree-Corrected Stochastic Blockmodel
#> -------------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 1
#> 
#> Traditional DCSBM parameterization:
#> 
#> Block memberships (z): 100 [factor] 
#> Degree heterogeneity (theta): 100 [numeric] 
#> Block probabilities (pi): 1 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 1 [dgeMatrix] 
#> S: 1 x 1 [ddiMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 50
#> Expected degree: 0.5
#> Expected density: 0.01

theta <- round(stats::rlnorm(100, 2))

cl2 <- chung_lu(
  theta = theta,
  expected_degree = 5
)

cl2
#> Undirected Degree-Corrected Stochastic Blockmodel
#> -------------------------------------------------
#> 
#> Nodes (n): 100 (arranged by block)
#> Blocks (k): 1
#> 
#> Traditional DCSBM parameterization:
#> 
#> Block memberships (z): 100 [factor] 
#> Degree heterogeneity (theta): 100 [numeric] 
#> Block probabilities (pi): 1 [numeric] 
#> 
#> Factor model parameterization:
#> 
#> X: 100 x 1 [dgeMatrix] 
#> S: 1 x 1 [ddiMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 500
#> Expected degree: 5
#> Expected density: 0.10101

edgelist <- sample_edgelist(cl)
edgelist
#> # A tibble: 28 × 2
#>     from    to
#>    <int> <int>
#>  1    14    76
#>  2     3     4
#>  3     8    23
#>  4     5     8
#>  5    20    59
#>  6     1     2
#>  7     8    33
#>  8     1    80
#>  9     2    43
#> 10     1     5
#> # ℹ 18 more rows
```
