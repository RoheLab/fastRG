# Create a directed factor model graph

A directed factor model graph is a directed generalized Poisson random
dot product graph. The edges in this graph are assumpted to be
independent and Poisson distributed. The graph is parameterized by its
expected adjacency matrix, with is `E[A] = X S Y'`. We do not recommend
that causal users use this function, see instead
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md)
and related functions, which will formulate common variants of the
stochastic blockmodels as undirected factor models *with lots of helpful
input validation*.

## Usage

``` r
directed_factor_model(
  X,
  S,
  Y,
  ...,
  expected_in_degree = NULL,
  expected_out_degree = NULL,
  expected_density = NULL,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
)
```

## Arguments

- X:

  A [`matrix()`](https://rdrr.io/r/base/matrix.html) or
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  representing real-valued latent node positions encoding community
  structure of incoming edges. Entries must be positive.

- S:

  A [`matrix()`](https://rdrr.io/r/base/matrix.html) or
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  mixing matrix. Entries must be positive.

- Y:

  A [`matrix()`](https://rdrr.io/r/base/matrix.html) or
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)
  representing real-valued latent node positions encoding community
  structure of outgoing edges. Entries must be positive.

- ...:

  Ignored. For internal developer use only.

- expected_in_degree:

  If specified, the desired expected in degree of the graph. Specifying
  `expected_in_degree` simply rescales `S` to achieve this. Defaults to
  `NULL`. Specify only one of `expected_in_degree`,
  `expected_out_degree`, and `expected_density`.

- expected_out_degree:

  If specified, the desired expected out degree of the graph. Specifying
  `expected_out_degree` simply rescales `S` to achieve this. Defaults to
  `NULL`. Specify only one of `expected_in_degree`,
  `expected_out_degree`, and `expected_density`.

- expected_density:

  If specified, the desired expected density of the graph. Specifying
  `expected_density` simply rescales `S` to achieve this. Defaults to
  `NULL`. Specify only one of `expected_in_degree`,
  `expected_out_degree`, and `expected_density`.

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

## Examples

``` r
n <- 1000

k1 <- 5
k2 <- 3

d <- 500

X <- matrix(rpois(n = n * k1, 1), nrow = n)
S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
Y <- matrix(rexp(n = k2 * d, 1), nrow = d)

fm <- directed_factor_model(X, S, Y)
fm
#> Directed Factor Model
#> ---------------------
#> 
#> Incoming Nodes (n): 1000
#> Incoming Rank (k1): 5
#> Outgoing Rank (k2): 3
#> Outgoing Nodes (d): 500
#> 
#> X: 1000 x 5 [dgeMatrix] 
#> S: 5 x 3 [dgeMatrix] 
#> Y: 500 x 3 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 365525
#> Expected density: 0.73105
#> Expected in degree: 731.1
#> Expected out degree: 365.5

fm2 <- directed_factor_model(X, S, Y, expected_in_degree = 50)
fm2
#> Directed Factor Model
#> ---------------------
#> 
#> Incoming Nodes (n): 1000
#> Incoming Rank (k1): 5
#> Outgoing Rank (k2): 3
#> Outgoing Nodes (d): 500
#> 
#> X: 1000 x 5 [dgeMatrix] 
#> S: 5 x 3 [dgeMatrix] 
#> Y: 500 x 3 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 25000
#> Expected density: 0.05
#> Expected in degree: 50
#> Expected out degree: 25
```
