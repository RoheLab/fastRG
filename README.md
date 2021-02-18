
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastRG

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/RoheLab/fastRG/branch/master/graph/badge.svg)](https://codecov.io/gh/RoheLab/fastRG?branch=master)
[![R build
status](https://github.com/RoheLab/fastRG/workflows/R-CMD-check/badge.svg)](https://github.com/RoheLab/fastRG/actions)
<!-- badges: end -->

`fastRG` quickly samples a broad class of network models known as
generalized random dot product graphs (GRDPGs). In particular, for
matrices *X*, *S* and *Y*, `fastRG` samples a matrix *A* with
expectation *X**S**Y*<sup>*T*</sup> where the entries are independently
Poisson distributed conditional on *X* and *Y*. This is primarily useful
when *A* is the adjacency matrix of a graph. Crucially, the sampling is
𝒪(*m*), where *m* is the number of the edges in graph, as opposed to the
naive sampling approach, which is 𝒪(*n*<sup>2</sup>), where *n* is the
number of nodes in the network.

For additional details, see the
[paper](https://arxiv.org/abs/1703.02998).

## Installation

`fastRG` is not yet on CRAN. You can install the development version
with:

``` r
# install.package("devtools")
devtools::install_github("RoheLab/fastRG")
```

## Usage

There are two stages to sampling from generalized random dot product
graphs. First, we sample the latent factors *X* and *Y*. Then we sample
*A* conditional on those latent factors. `fastRG` mimics this two-stage
sample structure. For example, to sample from a stochastic blockmodel,
we first create the latent factors.

``` r
library(fastRG)
#> Loading required package: Matrix

set.seed(27)

sbm <- sbm(n = 1000, k = 5, expected_density = 0.01)
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.
```

You can specify the latent factors and the mixing matrix *B* yourself,
but there are also defaults to enable fast prototyping. Here *B* was
randomly generated with `Uniform[0, 1]` entries and nodes were assigned
randomly to communities with equal probability of falling in all
communities. Printing the result object gives us some additional
information:

``` r
sbm
#> Undirected Stochastic Blockmodel
#> --------------------------------
#> 
#> Nodes (n): 1000 (arranged by block)
#> Blocks (k): 5
#> 
#> Traditional SBM parameterization:
#> 
#> Block memberships (z): 1000 [factor] 
#> Block probabilities (pi): 5 [numeric] 
#> Edge distribution: poisson
#> 
#> Factor model parameterization:
#> 
#> X: 1000 x 5 [dgCMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Expected edges: 10000
#> Expected degree: 10
#> Expected density: 0.01
```

Now, conditional on this latent representation, we can sample graphs.
`fastRG` supports several different output types, each of which is
specified by the suffix to `sample_*()` functions. For example, we can
obtain an edgelist in a `tibble` with:

``` r
sample_edgelist(sbm)
#> # A tibble: 9,986 x 2
#>     from    to
#>    <int> <int>
#>  1     1   184
#>  2   111     4
#>  3    86    12
#>  4    43   194
#>  5    61    37
#>  6    22    16
#>  7     4    10
#>  8    30   107
#>  9   119   209
#> 10    41    91
#> # … with 9,976 more rows
```

but we can just as easily obtain the graph as a sparse matrix

``` r
A <- sample_sparse(sbm)
A[1:10, 1:10]
#> 10 x 10 sparse Matrix of class "dsCMatrix"
#>                          
#>  [1,] . . . 1 . . . . . .
#>  [2,] . . 1 . . . . . . .
#>  [3,] . 1 . . . . . . . .
#>  [4,] 1 . . . . . . . . .
#>  [5,] . . . . . . . . . .
#>  [6,] . . . . . . . . . .
#>  [7,] . . . . . . . . . .
#>  [8,] . . . . . . . . . .
#>  [9,] . . . . . . . . . .
#> [10,] . . . . . . . . . .
```

or an igraph object

``` r
sample_igraph(sbm)
#> IGRAPH 07e1a72 U--- 1000 10072 -- 
#> + edges from 07e1a72:
#>  [1]  50-- 54  90--202  94--165 115--210 171--189   7-- 38  86--184  91--215
#>  [9]   3-- 92  35--111   9--159   3--158  66--131  11-- 28 164--214  48--163
#> [17] 116--141   8--189  10--170 102--193  14--207  99--209  36-- 89  72--213
#> [25]  62--126  17--136  14--145  15-- 15  41--211 161--174  23--215   5--132
#> [33] 121--190  61-- 84  29-- 95  28--182  52--109  22-- 56  24--166  14--109
#> [41]   2-- 37  93--125  53--153  27-- 62 122--211  52-- 81  34--180  12-- 93
#> [49] 120--196  50--190   7--175 107--186 137--155  92--173 195--213   2--163
#> [57] 159--208  98--215  53--205  75--124 164--167  27-- 44  82-- 96 177--188
#> [65]  55--167  38--134 170--192   5--  6  14--145 166--209  40--183  66--138
#> + ... omitted several edges
```

Note that every time we call `sample_*()` we draw a new sample.

``` r
A <- sample_sparse(sbm)
B <- sample_sparse(sbm)

all(A == B)  # random realizations from the SBM don't match!
#> [1] FALSE
```

## Efficient spectral decompositions

If you would like to obtain the singular value decomposition of the
population adjacency matrix conditional on latent factors, that is
straightforward:

``` r
s <- eigs_sym(sbm)
s$values
#> [1] 10.210177  3.676750  1.337299 -1.049310 -1.623513
```

Note that eigendecompositions and SVDS (for directed graphs) use
`RSpectra` and do not require explicitly forming large dense population
adjacency matrices; the population decompositions should be efficient in
both time and space for even large graphs.

## Key sampling options

There are several essential tools to modify graph sampling that you
should know about. First there are options that affect the latent factor
sampling:

-   `expected_degree`: Set the expected average degree of the graph by
    scaling sampling probabilities. We *strongly, strongly* recommend
    that you always set this option. If you do not, it is easy
    accidentally sample from large and dense graphs.

-   `expected_density`: Set the expected density of the graph by scaling
    sampling probabilities. You cannot specify both `expected_degree`
    and `expected_density` at the same time.

In the second stage of graph sampling, the options are:

-   `poisson_edges`: Either `TRUE` or `FALSE` depending on whether you
    would like a Bernoulli graph or a Poisson multi-graph. Scaling via
    `expected_degree` assumes a Poisson multi-graph, with some limited
    exceptions.

-   `allow_self_edges`: Whether nodes should be allowed to connect to
    themselves. Either `TRUE` or `FALSE`.

## Related work

You can find the original research code associated with `fastRG`
\[here\]\[<https://github.com/raningtky/sampleRDPG>\]. There is also a
Python translation of the original code in Python
[here](https://github.com/yunjhongwu/matrix-routines/blob/master/fastRG.py).
Both of these implementations are bare bones.
