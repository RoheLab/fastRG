# Sample a random dot product graph as a tidygraph graph

There are two steps to using the `fastRG` package. First, you must
parameterize a random dot product graph by sampling the latent factors.
Use functions such as
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md), etc,
to perform this specification. Then, use `sample_*()` functions to
generate a random graph in your preferred format.

## Usage

``` r
sample_tidygraph(factor_model, ...)

# S3 method for class 'undirected_factor_model'
sample_tidygraph(factor_model, ...)

# S3 method for class 'directed_factor_model'
sample_tidygraph(factor_model, ...)
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

A
[`tidygraph::tbl_graph()`](https://tidygraph.data-imaginist.com/reference/tbl_graph.html)
object that is possibly a multigraph (that is, we take there to be
multiple edges rather than weighted edges).

When `factor_model` is **undirected**:

    - the graph is undirected and one-mode.

When `factor_model` is **directed** and **square**:

    - the graph is directed and one-mode.

When `factor_model` is **directed** and **rectangular**:

    - the graph is undirected and bipartite.

Note that working with bipartite graphs in `tidygraph` is more complex
than working with one-mode graphs.

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
[`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md),
[`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md)

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
#> IGRAPH 90a0433 UN-- 100 213 -- 
#> + attr: name (v/c)
#> + edges from 90a0433 (vertex names):
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
#> # A tibble: 95 × 2
#>     from    to
#>    <int> <int>
#>  1    99    16
#>  2    88    15
#>  3     8    37
#>  4    65    48
#>  5    25     9
#>  6    19    47
#>  7     9    49
#>  8    97    19
#>  9    75    50
#> 10    36    40
#> # ℹ 85 more rows

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
#> # A tbl_graph: 150 nodes and 88 edges
#> #
#> # A bipartite multigraph with 69 components
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
#> # Edge Data: 88 × 2
#>    from    to
#>   <int> <int>
#> 1    83   101
#> 2    99   101
#> 3    90   102
#> # ℹ 85 more rows
```
