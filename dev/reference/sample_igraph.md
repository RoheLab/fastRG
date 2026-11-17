# Sample a random dot product graph as an igraph graph

There are two steps to using the `fastRG` package. First, you must
parameterize a random dot product graph by sampling the latent factors.
Use functions such as
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md), etc,
to perform this specification. Then, use `sample_*()` functions to
generate a random graph in your preferred format.

## Usage

``` r
sample_igraph(factor_model, ..., multigraph = TRUE)

# S3 method for class 'undirected_factor_model'
sample_igraph(factor_model, ..., multigraph = TRUE)

# S3 method for class 'directed_factor_model'
sample_igraph(factor_model, ..., multigraph = TRUE)
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

An
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
object. For Poisson graphs (`poisson_edges = TRUE`) with
`multigraph = TRUE` (the default), this is possibly a multigraph (that
is, we take there to be multiple edges rather than weighted edges). With
`multigraph = FALSE`, repeated edges are instead collapsed into a single
edge, and every edge has an integer `weight` attribute equal to the
corresponding entry of the
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md)
adjacency matrix (all weights are `1` when `poisson_edges = FALSE`).

When `factor_model` is **undirected**:

    - the graph is undirected and one-mode.

When `factor_model` is **directed** and **square**:

    - the graph is directed and one-mode.

When `factor_model` is **directed** and **rectangular**:

    - the graph is undirected and bipartite.

Note that working with bipartite graphs in `igraph` is more complex than
working with one-mode graphs.

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
[`sample_edgelist()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.md),
[`sample_edgelist.matrix()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.matrix.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md),
[`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md)

## Examples

``` r

library(igraph)
library(tidygraph)

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
#> Expected edges: 506
#> Expected degree: 10
#> Expected density: 0.1

### sampling graphs as edgelists ----------------------

edgelist <- sample_edgelist(ufm)
edgelist
#> # A tibble: 508 × 2
#>     from    to
#>    <int> <int>
#>  1    66    71
#>  2    83    87
#>  3    37    87
#>  4    92    95
#>  5    14    94
#>  6    54    85
#>  7    11    76
#>  8    45    57
#>  9    57    74
#> 10    11    22
#> # ℹ 498 more rows

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
#> IGRAPH 7382a91 UN-- 100 528 -- 
#> + attr: name (v/c)
#> + edges from 7382a91 (vertex names):
#>  [1] 15--36  61--87  65--80  57--71  3 --70  16--85  6 --16  1 --64  22--23 
#> [10] 52--85  6 --50  30--59  92--94  54--65  24--85  55--66  83--87  14--80 
#> [19] 32--54  12--19  13--30  6 --44  55--79  24--70  22--73  13--61  46--66 
#> [28] 24--83  12--37  46--73  73--84  61--87  25--75  16--47  15--66  17--59 
#> [37] 25--75  17--74  1 --85  67--87  25--80  23--88  13--100 14--27  1 --84 
#> [46] 43--85  47--79  61--82  22--63  26--76  14--16  22--29  4 --25  78--92 
#> [55] 69--87  84--84  55--87  19--95  3 --82  86--94  98--99  10--22  66--70 
#> [64] 38--94  79--91  12--25  45--79  47--56  18--64  13--78  22--90  37--62 
#> + ... omitted several edges

### sampling graphs as tidygraph graphs ---------------

sample_tidygraph(ufm)
#> # A tbl_graph: 100 nodes and 508 edges
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
#> # Edge Data: 508 × 2
#>    from    to
#>   <int> <int>
#> 1    17    44
#> 2    47    94
#> 3    19    80
#> # ℹ 505 more rows

### weighted graphs instead of multigraphs ------------

# collapse repeated edges into a single edge with an integer weight

sample_edgelist(ufm, multigraph = FALSE)
#> # A tibble: 478 × 3
#>     from    to weight
#>    <int> <int>  <int>
#>  1     1    23      2
#>  2     1    47      1
#>  3     1    49      1
#>  4     1    54      1
#>  5     1    65      1
#>  6     1    73      1
#>  7     1    83      1
#>  8     1    84      1
#>  9     1    85      1
#> 10     1    99      1
#> # ℹ 468 more rows
sample_igraph(ufm, multigraph = FALSE)
#> IGRAPH ee61bd5 UNW- 100 502 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from ee61bd5 (vertex names):
#>  [1] 1--2  1--11 1--27 1--32 1--52 1--55 1--58 1--62 1--70 1--85 1--99 2--4 
#> [13] 2--9  2--13 2--17 2--25 2--30 2--40 2--47 2--56 2--68 2--71 2--79 2--93
#> [25] 3--13 3--29 3--44 3--45 3--52 3--58 3--60 3--68 3--71 3--81 3--90 3--92
#> [37] 3--99 4--6  4--16 4--37 4--68 4--74 4--77 4--88 4--97 5--9  5--11 5--12
#> [49] 5--21 5--25 5--29 5--38 5--41 5--44 5--57 5--69 5--73 5--75 5--79 5--82
#> [61] 5--84 5--91 6--36 6--66 6--73 7--11 7--16 7--22 7--66 8--47 8--53 8--64
#> [73] 8--66 8--71 8--85 9--10 9--25 9--35 9--45 9--50 9--55 9--58 9--61 9--62
#> [85] 9--65 9--73 9--74 9--75 9--77 9--79 9--85 9--91 9--94 9--97
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
#>  1    72    23
#>  2    47    15
#>  3    27    32
#>  4    85     8
#>  5    87    23
#>  6    52    49
#>  7    77    45
#>  8    95    10
#>  9    38    34
#> 10    56     8
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
#> # A tbl_graph: 150 nodes and 109 edges
#> #
#> # A bipartite multigraph with 56 components
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
#> # Edge Data: 109 × 2
#>    from    to
#>   <int> <int>
#> 1    34   101
#> 2    37   101
#> 3    46   101
#> # ℹ 106 more rows
```
