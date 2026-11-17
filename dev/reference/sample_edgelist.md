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
sample_edgelist(factor_model, ..., multigraph = TRUE)

# S3 method for class 'undirected_factor_model'
sample_edgelist(factor_model, ..., multigraph = TRUE)

# S3 method for class 'directed_factor_model'
sample_edgelist(factor_model, ..., multigraph = TRUE)
```

## Arguments

- factor_model:

  A
  [`directed_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/directed_factor_model.md)
  or
  [`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md).

- ...:

  Ignored. Do not use.

- multigraph:

  Logical scalar, only relevant for Poisson graphs
  (`poisson_edges = TRUE`), where a pair of nodes can be connected by
  more than one edge. When `TRUE` (the default), such repeated edges are
  represented as a multigraph, that is, as repeated edges (or repeated
  rows of an edgelist) without any edge weights. When `FALSE`, repeated
  edges between the same pair of nodes are collapsed into a single edge
  with an integer `weight` attribute that counts the number of times the
  edge was sampled. This weight is the corresponding entry of the
  adjacency matrix returned by
  [`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md).
  In undirected graphs, the edges `(i, j)` and `(j, i)` are the same
  edge and are collapsed together. When `multigraph = FALSE` the
  `weight` attribute is always present, even for Bernoulli graphs
  (`poisson_edges = FALSE`), where every weight is `1`.

## Value

A single realization of a random Poisson (or Bernoulli) Dot Product
Graph, represented as a
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with two integer columns, `from` and `to`. When `multigraph = FALSE`,
there is an additional integer column `weight` and each `(from, to)`
pair appears at most once.

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
#> # A tibble: 256 × 2
#>     from    to
#>    <int> <int>
#>  1    44    57
#>  2    65    66
#>  3    75    97
#>  4    12    87
#>  5    13    71
#>  6    24    87
#>  7    37    65
#>  8    92    98
#>  9    14    94
#> 10    22    85
#> # ℹ 246 more rows

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
#> IGRAPH 7256255 UN-- 100 213 -- 
#> + attr: name (v/c)
#> + edges from 7256255 (vertex names):
#>  [1] 40--56  71--85  80--80  25--97  66--94  6 --22  47--86  15--50  62--94 
#> [10] 56--98  12--75  26--30  14--84  13--95  85--94  15--40  47--95  17--70 
#> [19] 46--78  44--84  6 --68  94--98  26--38  61--87  48--100 6 --36  12--93 
#> [28] 31--85  56--97  62--83  1 --55  3 --75  24--32  16--85  25--85  24--71 
#> [37] 9 --65  2 --71  25--85  46--100 44--87  9 --83  1 --46  44--61  84--95 
#> [46] 84--95  5 --12  15--23  64--73  55--58  25--59  12--25  68--85  79--85 
#> [55] 61--75  46--62  13--95  50--58  17--87  11--13  85--100 10--14  27--74 
#> [64] 20--64  3 --15  47--52  11--58  86--87  3 --43  20--59  27--99  35--86 
#> + ... omitted several edges

### sampling graphs as tidygraph graphs ---------------

sample_tidygraph(ufm)
#> # A tbl_graph: 100 nodes and 239 edges
#> #
#> # An undirected multigraph with 2 components
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
#> # Edge Data: 239 × 2
#>    from    to
#>   <int> <int>
#> 1    16    66
#> 2    77    87
#> 3    27    76
#> # ℹ 236 more rows

### weighted graphs instead of multigraphs ------------

# collapse repeated edges into a single edge with an integer weight

sample_edgelist(ufm, multigraph = FALSE)
#> # A tibble: 263 × 3
#>     from    to weight
#>    <int> <int>  <int>
#>  1     1    37      1
#>  2     1    49      1
#>  3     1    74      1
#>  4     2     5      1
#>  5     2    30      1
#>  6     2    38      1
#>  7     2    61      1
#>  8     2    64      1
#>  9     2    72      1
#> 10     2    74      1
#> # ℹ 253 more rows
sample_igraph(ufm, multigraph = FALSE)
#> IGRAPH 5391c07 UNW- 100 245 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from 5391c07 (vertex names):
#>  [1] 1 --17  1 --56  1 --99  2 --20  2 --24  2 --63  2 --70  2 --87  2 --100
#> [10] 3 --59  3 --73  3 --87  4 --7   4 --8   4 --31  5 --16  5 --80  5 --100
#> [19] 6 --61  6 --63  6 --70  6 --98  8 --14  9 --14  9 --19  9 --24  9 --37 
#> [28] 9 --38  9 --44  9 --46  9 --59  9 --76  9 --77  9 --78  9 --85  9 --94 
#> [37] 10--10  10--20  10--31  10--49  10--69  10--75  10--86  11--58  11--59 
#> [46] 11--63  12--19  12--49  12--54  12--69  12--75  12--82  13--45  13--47 
#> [55] 13--76  13--95  14--59  14--69  15--34  16--16  16--46  16--63  16--84 
#> [64] 16--96  16--98  17--29  17--31  17--36  17--46  17--49  17--52  17--62 
#> + ... omitted several edges

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
#> # A tibble: 109 × 2
#>     from    to
#>    <int> <int>
#>  1    94    38
#>  2    95     5
#>  3    12    10
#>  4    61     6
#>  5    40    17
#>  6    47    46
#>  7    54    27
#>  8    24     5
#>  9    19    41
#> 10    40     2
#> # ℹ 99 more rows

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
#> # A tbl_graph: 150 nodes and 98 edges
#> #
#> # A bipartite multigraph with 65 components
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
#> # Edge Data: 98 × 2
#>    from    to
#>   <int> <int>
#> 1    75   102
#> 2    10   104
#> 3    79   104
#> # ℹ 95 more rows
```
