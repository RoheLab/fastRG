# Sample a random edgelist from a random dot product graph

There are two steps to using the `fastRG` package. First, you must
parameterize a random dot product graph by sampling the latent factors.
Use functions such as
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md), etc,
to perform this specification. Then, use `sample_*()` functions to
generate a random graph in your preferred format.

## Usage

``` r
sample_edgelist(factor_model, ...)

# S3 method for class 'undirected_factor_model'
sample_edgelist(factor_model, ...)

# S3 method for class 'directed_factor_model'
sample_edgelist(factor_model, ...)
```

## Arguments

- factor_model:

  A
  [`directed_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/directed_factor_model.md)
  or
  [`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md).

- ...:

  Ignored. Do not use.

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
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md),
and
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md)
over `sample_edgelist()`.

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
[`sample_edgelist.matrix()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.matrix.md),
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md)

## Examples

``` r
library(igraph)
#> 
#> Attaching package: ‘igraph’
#> The following objects are masked from ‘package:fastRG’:
#> 
#>     chung_lu, sbm
#> The following objects are masked from ‘package:stats’:
#> 
#>     decompose, spectrum
#> The following object is masked from ‘package:base’:
#> 
#>     union
library(tidygraph)
#> 
#> Attaching package: ‘tidygraph’
#> The following object is masked from ‘package:igraph’:
#> 
#>     groups
#> The following object is masked from ‘package:stats’:
#> 
#>     filter

set.seed(27)

##### undirected examples ----------------------------

n <- 100
k <- 5

X <- matrix(rpois(n = n * k, 1), nrow = n)
S <- matrix(runif(n = k * k, 0, .1), nrow = k)

# S will be symmetrized internal here, or left unchanged if
# it is already symmetric

ufm <- undirected_factor_model(
  X, S,
  expected_density = 0.1
)

ufm
#> Undirected Factor Model
#> -----------------------
#> 
#> Nodes (n): 100
#> Rank (k): 5
#> 
#> X: 100 x 5 [dgeMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 495
#> Expected degree: 5
#> Expected density: 0.1

### sampling graphs as edgelists ----------------------

edgelist <- sample_edgelist(ufm)
edgelist
#> # A tibble: 500 × 2
#>     from    to
#>    <int> <int>
#>  1    66    71
#>  2    85    87
#>  3    37    54
#>  4    70    92
#>  5    14    44
#>  6    66    85
#>  7    76    83
#>  8    57    87
#>  9    57    95
#> 10    22    94
#> # ℹ 490 more rows

### sampling graphs as sparse matrices ----------------

A <- sample_sparse(ufm)

inherits(A, "dsCMatrix")
#> [1] TRUE
isSymmetric(A)
#> [1] TRUE
dim(A)
#> [1] 100 100

B <- sample_sparse(ufm)

inherits(B, "dsCMatrix")
#> [1] TRUE
isSymmetric(B)
#> [1] TRUE
dim(B)
#> [1] 100 100

### sampling graphs as igraph graphs ------------------

sample_igraph(ufm)
#> IGRAPH 0199eee UN-- 100 486 -- 
#> + attr: name (v/c)
#> + edges from 0199eee (vertex names):
#>  [1] 65--87  84--100 12--87  13--95  3 --92  25--94  54--98  16--22  1 --66 
#> [10] 13--94  65--79  12--66  79--94  55--56  30--64  13--22  22--40  37--80 
#> [19] 88--95  11--22  85--94  52--94  11--37  12--16  19--75  47--74  63--97 
#> [28] 12--61  11--73  2 --71  25--28  61--70  88--98  44--71  61--97  46--56 
#> [37] 14--85  36--65  14--17  20--71  12--12  57--85  59--71  46--90  30--38 
#> [46] 17--55  59--98  15--47  37--62  49--85  65--98  37--98  22--33  56--77 
#> [55] 25--51  20--80  16--57  25--71  52--64  12--47  8 --80  18--79  22--62 
#> [64] 14--31  37--69  16--54  26--90  38--94  20--79  70--97  19--90  11--71 
#> + ... omitted several edges

### sampling graphs as tidygraph graphs ---------------

sample_tidygraph(ufm)
#> # A tbl_graph: 100 nodes and 501 edges
#> #
#> # An undirected multigraph with 1 component
#> #
#> # Node Data: 100 × 1 (active)
#>    name 
#>    <chr>
#>  1 1    
#>  2 2    
#>  3 3    
#>  4 4    
#>  5 5    
#>  6 6    
#>  7 7    
#>  8 8    
#>  9 9    
#> 10 10   
#> # ℹ 90 more rows
#> #
#> # Edge Data: 501 × 2
#>    from    to
#>   <int> <int>
#> 1    54    94
#> 2    56    94
#> 3    16    22
#> # ℹ 498 more rows

##### directed examples ----------------------------

n2 <- 100

k1 <- 5
k2 <- 3

d <- 50

X <- matrix(rpois(n = n2 * k1, 1), nrow = n2)
S <- matrix(runif(n = k1 * k2, 0, .1), nrow = k1, ncol = k2)
Y <- matrix(rexp(n = k2 * d, 1), nrow = d)

fm <- directed_factor_model(X, S, Y, expected_in_degree = 2)
fm
#> Directed Factor Model
#> ---------------------
#> 
#> Incoming Nodes (n): 100
#> Incoming Rank (k1): 5
#> Outgoing Rank (k2): 3
#> Outgoing Nodes (d): 50
#> 
#> X: 100 x 5 [dgeMatrix] 
#> S: 5 x 3 [dgeMatrix] 
#> Y: 50 x 3 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 100
#> Expected density: 0.02
#> Expected in degree: 2
#> Expected out degree: 1

### sampling graphs as edgelists ----------------------

edgelist2 <- sample_edgelist(fm)
edgelist2
#> # A tibble: 105 × 2
#>     from    to
#>    <int> <int>
#>  1    84    34
#>  2    80    16
#>  3    42    30
#>  4    42    31
#>  5    47    26
#>  6     7    31
#>  7    39    14
#>  8    14    11
#>  9    49    47
#> 10    54    28
#> # ℹ 95 more rows

### sampling graphs as sparse matrices ----------------

A2 <- sample_sparse(fm)

inherits(A2, "dgCMatrix")
#> [1] TRUE
isSymmetric(A2)
#> [1] FALSE
dim(A2)
#> [1] 100  50

B2 <- sample_sparse(fm)

inherits(B2, "dgCMatrix")
#> [1] TRUE
isSymmetric(B2)
#> [1] FALSE
dim(B2)
#> [1] 100  50

### sampling graphs as igraph graphs ------------------

# since the number of rows and the number of columns
# in `fm` differ, we will get a bipartite igraph here

# creating the bipartite igraph is slow relative to other
# sampling -- if this is a blocker for
# you please open an issue and we can investigate speedups

dig <- sample_igraph(fm)
is_bipartite(dig)
#> [1] TRUE

### sampling graphs as tidygraph graphs ---------------

sample_tidygraph(fm)
#> # A tbl_graph: 150 nodes and 105 edges
#> #
#> # A bipartite multigraph with 59 components
#> #
#> # Node Data: 150 × 1 (active)
#>    type 
#>    <lgl>
#>  1 FALSE
#>  2 FALSE
#>  3 FALSE
#>  4 FALSE
#>  5 FALSE
#>  6 FALSE
#>  7 FALSE
#>  8 FALSE
#>  9 FALSE
#> 10 FALSE
#> # ℹ 140 more rows
#> #
#> # Edge Data: 105 × 2
#>    from    to
#>   <int> <int>
#> 1    55   101
#> 2    61   101
#> 3    89   101
#> # ℹ 102 more rows
```
