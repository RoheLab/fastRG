# Create an directed erdos renyi object

Create an directed erdos renyi object

## Usage

``` r
directed_erdos_renyi(
  n,
  ...,
  p = NULL,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
)
```

## Arguments

- n:

  Number of nodes in graph.

- ...:

  Arguments passed on to
  [`directed_factor_model`](https://rohelab.github.io/fastRG/reference/directed_factor_model.md)

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

- p:

  Probability of an edge between any two nodes. You must specify either
  `p`, `expected_in_degree`, or `expected_out_degree`.

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

A `directed_factor_model` S3 class based on a list with the following
elements:

- `X`: The incoming latent positions as a
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  object.

- `S`: The mixing matrix as a
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  object.

- `Y`: The outgoing latent positions as a
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  object.

- `n`: The number of nodes with incoming edges in the network.

- `k1`: The dimension of the latent node position vectors encoding
  incoming latent communities (i.e. in `X`).

- `d`: The number of nodes with outgoing edges in the network. Does not
  need to match `n` – rectangular adjacency matrices are supported.

- `k2`: The dimension of the latent node position vectors encoding
  outgoing latent communities (i.e. in `Y`).

- `poisson_edges`: Whether or not the graph is taken to be have Poisson
  or Bernoulli edges, as indicated by a logical vector of length 1.

- `allow_self_loops`: Whether or not self loops are allowed.

## See also

Other erdos renyi:
[`erdos_renyi()`](https://rohelab.github.io/fastRG/reference/erdos_renyi.md)

Other directed graphs:
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/reference/directed_dcsbm.md)

## Examples

``` r
set.seed(87)

er <- directed_erdos_renyi(n = 10, p = 0.1)
er
#> Directed Factor Model
#> ---------------------
#> 
#> Incoming Nodes (n): 10
#> Incoming Rank (k1): 1
#> Outgoing Rank (k2): 1
#> Outgoing Nodes (d): 10
#> 
#> X: 10 x 1 [dgeMatrix] 
#> S: 1 x 1 [ddiMatrix] 
#> Y: 10 x 1 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 10
#> Expected density: 0.1
#> Expected in degree: 1
#> Expected out degree: 1

big_er <- directed_erdos_renyi(n = 1000, expected_in_degree = 5)
big_er
#> Directed Factor Model
#> ---------------------
#> 
#> Incoming Nodes (n): 1000
#> Incoming Rank (k1): 1
#> Outgoing Rank (k2): 1
#> Outgoing Nodes (d): 1000
#> 
#> X: 1000 x 1 [dgeMatrix] 
#> S: 1 x 1 [ddiMatrix] 
#> Y: 1000 x 1 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 5000
#> Expected density: 0.005
#> Expected in degree: 5
#> Expected out degree: 5

A <- sample_sparse(er)
A
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                          
#>  [1,] . . . . . . . . . .
#>  [2,] . . . . . . 1 . . .
#>  [3,] . . . 1 . . . . . .
#>  [4,] . . . . 1 . 1 . 1 1
#>  [5,] . . . . . . . . . .
#>  [6,] . . . . . . 1 . . .
#>  [7,] . . . . . . . . . .
#>  [8,] 1 . . . . . . . . .
#>  [9,] . . . . . 1 . . . .
#> [10,] . . . . . . 2 . . .
```
