# Low level interface to sample RPDG edgelists

**This is a brakes-off, no safety checks interface.** We strongly
recommend that you do not call `sample_edgelist.matrix()` unless you
know what you are doing, and even then, we still do not recommend it, as
you will bypass all typical input validation. **extremely loud
coughing** All those who bypass input validation suffer foolishly at
their own hand. **extremely loud coughing**

## Usage

``` r
# S3 method for class 'matrix'
sample_edgelist(
  factor_model,
  S,
  Y,
  directed,
  poisson_edges,
  allow_self_loops,
  ...
)

# S3 method for class 'Matrix'
sample_edgelist(
  factor_model,
  S,
  Y,
  directed,
  poisson_edges,
  allow_self_loops,
  ...
)
```

## Arguments

- factor_model:

  An `n` by `k1` [`matrix()`](https://rdrr.io/r/base/matrix.html) or
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html) of
  latent node positions encoding incoming edge community membership. The
  `X` matrix in Rohe et al (2017). Naming differs only for consistency
  with the S3 generic.

- S:

  A `k1` by `k2` mixing [`matrix()`](https://rdrr.io/r/base/matrix.html)
  or [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html).
  In the undirect case this is assumed to be symmetric but **we do not
  check that this is the case**.

- Y:

  A `d` by `k2` [`matrix()`](https://rdrr.io/r/base/matrix.html) or
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html) of
  latent node positions encoding outgoing edge community membership.

- directed:

  Logical indicating whether or not the graph should be directed. When
  `directed = FALSE`, symmetrizes `S` internally. `Y = X` together with
  a symmetric `S` implies a symmetric expectation (although not
  necessarily an undirected graph). When `directed = FALSE`, samples a
  directed graph with symmetric expectation, and then adds edges until
  symmetry is achieved.

- poisson_edges:

  Whether or not to remove duplicate edges after sampling. See Section
  2.3 of Rohe et al. (2017) for some additional details. Defaults to
  `TRUE`.

- allow_self_loops:

  Logical indicating whether or not nodes should be allowed to form
  edges with themselves. Defaults to `TRUE`. When `FALSE`, sampling
  proceeds allowing self-loops, and these are then removed after the
  fact.

- ...:

  Ignored, for generic consistency only.

## Value

A single realization of a random Poisson (or Bernoulli) Dot Product
Graph, represented as a
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with two integer columns, `from` and `to`.

**NOTE**: Indices for isolated nodes will not appear in the edgelist!
This can lead to issues if you construct network objects from the
edgelist directly.

In the undirected case, `from` and `to` do not encode information about
edge direction, but we will always have `from <= to` for convenience of
edge identification.

To avoid handling such considerations yourself, we recommend using
[`sample_sparse()`](https://rohelab.github.io/fastRG/reference/sample_sparse.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/reference/sample_igraph.md),
and
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/reference/sample_tidygraph.md)
over
[`sample_edgelist()`](https://rohelab.github.io/fastRG/reference/sample_edgelist.md).

## Details

This function implements the `fastRG` algorithm as described in Rohe et
al (2017). Please see the paper (which is short and open access!!) for
details.

## References

Rohe, Karl, Jun Tao, Xintian Han, and Norbert Binkiewicz. 2017. "A Note
on Quickly Sampling a Sparse Matrix with Low Rank Expectation." Journal
of Machine Learning Research; 19(77):1-13, 2018.
<https://www.jmlr.org/papers/v19/17-128.html>

## See also

Other samplers:
[`sample_edgelist()`](https://rohelab.github.io/fastRG/reference/sample_edgelist.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/reference/sample_igraph.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/reference/sample_sparse.md),
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/reference/sample_tidygraph.md)

## Examples

``` r
set.seed(46)

n <- 10000
d <- 1000

k1 <- 5
k2 <- 3

X <- matrix(rpois(n = n * k1, 1), nrow = n)
S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1)
Y <- matrix(rpois(n = d * k2, 1), nrow = d)

sample_edgelist(X, S, Y, TRUE, TRUE, TRUE)
#> # A tibble: 8,089,734 × 2
#>     from    to
#>    <int> <int>
#>  1  2325   310
#>  2   531   962
#>  3    83   779
#>  4  2894   294
#>  5  5454   157
#>  6  7379   677
#>  7  9789   178
#>  8  3077   750
#>  9  7587   992
#> 10  1929   387
#> # ℹ 8,089,724 more rows
```
