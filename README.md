
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastRG

<!-- badges: start -->

[![R-CMD-check](https://github.com/alexpghayes/fastRG/workflows/R-CMD-check/badge.svg)](https://github.com/alexpghayes/fastRG/actions)
[![Codecov test
coverage](https://codecov.io/gh/RoheLab/fastRG/branch/main/graph/badge.svg)](https://codecov.io/gh/RoheLab/fastRG?branch=main)
<!-- badges: end -->

`fastRG` quickly samples a broad class of network models known as
generalized random dot product graphs (GRDPGs). In particular, for
matrices *X*, *S* and *Y*, `fastRG` samples a matrix *A* with
expectation *X**S**Y*<sup>*T*</sup> where the entries are independently
Poisson distributed conditional on *X* and *Y*. This is primarily useful
when *A* is the adjacency matrix of a graph. Crucially, the sampling is
𝒪(*m*), where *m* is the number of the edges in graph, as opposed to the
naive sampling approach, which is 𝒪(*n*<sup>2</sup>), where *n* is the
number of nodes in the network. For additional details, see the
[paper](https://arxiv.org/abs/1703.02998).

`fastRG` has two primary use cases:

1.  Sampling enormous sparse graphs that cannot feasibly be sampled with
    existing samplers, and
2.  validating new methods for random dot product graphs (and variants).

`fastRG` makes the latent parameters of random dot product graphs
readily available to users, such that simulation studies for community
detection, subspace recovery, etc, are straightforward.

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
#> # A tibble: 4,990 x 2
#>     from    to
#>    <int> <int>
#>  1    94   111
#>  2    86   143
#>  3    43    89
#>  4    61   159
#>  5    22   210
#>  6     4   197
#>  7    30   145
#>  8   119   136
#>  9    41   142
#> 10   175   182
#> # … with 4,980 more rows
```

but we can just as easily obtain the graph as a sparse matrix

``` r
A <- sample_sparse(sbm)
A[1:10, 1:10]
#> 10 x 10 sparse Matrix of class "dsCMatrix"
#>                          
#>  [1,] . . . . . . . . . .
#>  [2,] . . . . . . . . . .
#>  [3,] . . . . . . . . . .
#>  [4,] . . . . . . . . . .
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
#> IGRAPH d4b2e7c UN-- 1000 5085 -- 
#> + attr: name (v/c)
#> + edges from d4b2e7c (vertex names):
#>  [1] 142--155 115--118 9  --82  20 --39  80 --127 134--216 62 --196 77 --199
#>  [9] 53 --209 1  --57  40 --98  1  --47  167--218 6  --188 172--179 88 --143
#> [17] 88 --185 28 --30  98 --124 39 --116 92 --118 160--183 35 --41  156--161
#> [25] 47 --104 92 --137 101--192 30 --125 45 --208 15 --109 135--203 86 --125
#> [33] 63 --166 4  --212 148--151 115--102 9  --92  61 --112 134--152 17 --83 
#> [41] 32 --146 41 --175 142--93  77 --208 8  --112 17 --171 83 --171 90 --204
#> [49] 86 --84  49 --140 22 --79  3  --205 5  --103 75 --133 78 --152 23 --70 
#> [57] 143--162 41 --12  49 --56  35 --99  78 --78  61 --170 9  --58  108--141
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
#> [1]  5.1050886  1.8383749  0.6686493 -0.5246550 -0.8117567
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

## Known issues

Sampling blockmodels with very small numbers of nodes (or blockmodels
with the number of blocks `k` on the same order as `n`) results in a
degeneracy that can cause issues.

## Related work

[`igraph`](https://igraph.org/r/) allows users to sample SBMs (in
𝒪(*m* + *n* + *k*<sup>2</sup>) time) and random dot product graphs (in
𝒪(*n*<sup>2</sup>*k*) time). You can find the original research code
associated with `fastRG`
[here](https://github.com/raningtky/sampleRDPG). There is also a Python
translation of the original code in Python
[here](https://github.com/yunjhongwu/matrix-routines/blob/master/fastRG.py).
Both of these implementations are bare bones.
