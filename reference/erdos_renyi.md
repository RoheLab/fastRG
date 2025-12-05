# Create an undirected erdos renyi object

Create an undirected erdos renyi object

## Usage

``` r
erdos_renyi(n, ..., p = NULL, poisson_edges = TRUE, allow_self_loops = TRUE)
```

## Arguments

- n:

  Number of nodes in graph.

- ...:

  Arguments passed on to
  [`undirected_factor_model`](https://rohelab.github.io/fastRG/reference/undirected_factor_model.md)

  `expected_degree`

  :   If specified, the desired expected degree of the graph. Specifying
      `expected_degree` simply rescales `S` to achieve this. Defaults to
      `NULL`. Do not specify both `expected_degree` and
      `expected_density` at the same time.

- p:

  Probability of an edge between any two nodes. You must specify either
  `p` or `expected_degree`.

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

An `undirected_factor_model` S3 class based on a list with the following
elements:

- `X`: The latent positions as a
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  object.

- `S`: The mixing matrix as a
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  object.

- `n`: The number of nodes in the network.

- `k`: The rank of expectation matrix. Equivalently, the dimension of
  the latent node position vectors.

## See also

Other erdos renyi:
[`directed_erdos_renyi()`](https://rohelab.github.io/fastRG/reference/directed_erdos_renyi.md)

Other undirected graphs:
[`chung_lu()`](https://rohelab.github.io/fastRG/reference/chung_lu.md),
[`dcsbm()`](https://rohelab.github.io/fastRG/reference/dcsbm.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/reference/planted_partition.md),
[`sbm()`](https://rohelab.github.io/fastRG/reference/sbm.md)

## Examples

``` r
set.seed(87)

er <- erdos_renyi(n = 10, p = 0.1)
er
#> Undirected Factor Model
#> -----------------------
#> 
#> Nodes (n): 10
#> Rank (k): 1
#> 
#> X: 10 x 1 [dgeMatrix] 
#> S: 1 x 1 [ddiMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 10
#> Expected degree: 1
#> Expected density: 0.22222


er <- erdos_renyi(n = 10, expected_density = 0.1)
er
#> Undirected Factor Model
#> -----------------------
#> 
#> Nodes (n): 10
#> Rank (k): 1
#> 
#> X: 10 x 1 [dgeMatrix] 
#> S: 1 x 1 [ddiMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 4
#> Expected degree: 0.4
#> Expected density: 0.1

big_er <- erdos_renyi(n = 1000, expected_degree = 5)
big_er
#> Undirected Factor Model
#> -----------------------
#> 
#> Nodes (n): 1000
#> Rank (k): 1
#> 
#> X: 1000 x 1 [dgeMatrix] 
#> S: 1 x 1 [ddiMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 5000
#> Expected degree: 5
#> Expected density: 0.01001

A <- sample_sparse(er)
A
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                          
#>  [1,] 0 . . . . . . . . .
#>  [2,] . 0 . . . . . . . .
#>  [3,] . . 0 . . . . . . .
#>  [4,] . . . 0 . . . . . .
#>  [5,] . . . . 0 . . . . .
#>  [6,] . . . . . 0 . . . .
#>  [7,] . . . . . . 0 . . .
#>  [8,] . . . . . . . 0 . .
#>  [9,] . . . . . . . . 0 .
#> [10,] . . . . . . . . . 0
```
