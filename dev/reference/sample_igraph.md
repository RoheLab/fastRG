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
sample_igraph(factor_model, ...)

# S3 method for class 'undirected_factor_model'
sample_igraph(factor_model, ...)

# S3 method for class 'directed_factor_model'
sample_igraph(factor_model, ...)
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

An
[`igraph::igraph()`](https://r.igraph.org/reference/aaa-igraph-package.html)
object that is possibly a multigraph (that is, we take there to be
multiple edges rather than weighted edges).

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
#> Expected edges: 495
#> Expected degree: 5
#> Expected density: 0.1

### sampling graphs as edgelists ----------------------

edgelist <- sample_edgelist(ufm)
edgelist
#> # A tibble: 251 × 2
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
#> # ℹ 241 more rows

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
#> IGRAPH 35ac6a0 UN-- 100 236 -- 
#> + attr: name (v/c)
#> + edges from 35ac6a0 (vertex names):
#>  [1] 16--98  45--98  57--88  19--88  22--100 1 --70  37--87  13--66  19--84 
#> [10] 17--74  37--92  65--83  1 --46  17--94  86--92  59--94  70--90  30--66 
#> [19] 10--61  44--87  43--79  39--92  46--73  6 --100 12--17  30--44  16--64 
#> [28] 17--54  11--71  64--100 9 --56  41--80  47--79  11--100 62--99  9 --14 
#> [37] 71--81  11--80  5 --11  32--76  1 --94  55--81  22--38  62--79  46--88 
#> [46] 41--82  37--79  76--98  73--86  17--75  34--74  52--56  71--90  2 --39 
#> [55] 14--14  17--61  3 --11  82--91  72--87  17--23  16--46  87--97  47--98 
#> [64] 61--66  5 --94  45--51  43--71  30--97  36--72  28--47  85--97  15--27 
#> + ... omitted several edges

### sampling graphs as tidygraph graphs ---------------

sample_tidygraph(ufm)
#> # A tbl_graph: 100 nodes and 238 edges
#> #
#> # An undirected multigraph with 8 components
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
#> # Edge Data: 238 × 2
#>    from    to
#>   <int> <int>
#> 1    46    80
#> 2    12    25
#> 3    54    66
#> # ℹ 235 more rows

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
#> # A tibble: 96 × 2
#>     from    to
#>    <int> <int>
#>  1    40    10
#>  2    46    36
#>  3    31    36
#>  4     4    34
#>  5    40     9
#>  6    54    32
#>  7    38    24
#>  8   100    17
#>  9    79    11
#> 10    53    11
#> # ℹ 86 more rows

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
#> # A tbl_graph: 150 nodes and 104 edges
#> #
#> # A bipartite multigraph with 58 components
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
#> # Edge Data: 104 × 2
#>    from    to
#>   <int> <int>
#> 1     2   102
#> 2    34   102
#> 3    45   102
#> # ℹ 101 more rows
```
