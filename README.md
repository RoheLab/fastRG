
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastRG

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/RoheLab/fastRG/branch/main/graph/badge.svg)](https://app.codecov.io/gh/RoheLab/fastRG?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/fastRG)](https://CRAN.R-project.org/package=fastRG)
[![R-CMD-check](https://github.com/RoheLab/fastRG/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RoheLab/fastRG/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`fastRG` quickly samples a broad class of network models known as
generalized random dot product graphs (GRDPGs). In particular, for
matrices $X$, $S$ and $Y$, `fastRG` samples a matrix $A$ with
expectation $X S Y^T$ where the entries are independently Poisson
distributed conditional on $X$ and $Y$. This is primarily useful when
$A$ is the adjacency matrix of a graph. Crucially, the sampling is
$\mathcal O(m)$, where $m$ is the number of the edges in graph, as
opposed to the naive sampling approach, which is $\mathcal O(n^2)$,
where $n$ is the number of nodes in the network. For additional details,
see the [paper](https://arxiv.org/abs/1703.02998) \[1\].

`fastRG` has two primary use cases:

1.  Sampling enormous sparse graphs that cannot feasibly be sampled with
    existing samplers, and
2.  validating new methods for random dot product graphs (and variants).

`fastRG` makes the latent parameters of random dot product graphs
readily available to users, such that simulation studies for community
detection, subspace recovery, etc, are straightforward.

## Installation

You can install the released version of fastRG from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fastRG")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("RoheLab/fastRG")
```

## Usage

There are two stages to sampling from generalized random dot product
graphs. First, we sample the latent factors $X$ and $Y$. Then we sample
$A$ conditional on those latent factors. `fastRG` mimics this two-stage
sample structure. For example, to sample from a stochastic blockmodel,
we first create the latent factors.

``` r
library(fastRG)
#> Loading required package: Matrix

set.seed(27)

sbm <- sbm(n = 1000, k = 5, expected_density = 0.01)
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.
```

You can specify the latent factors and the mixing matrix $B$ yourself,
but there are also defaults to enable fast prototyping. Here $B$ was
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
#> Factor model parameterization:
#> 
#> X: 1000 x 5 [dgCMatrix] 
#> S: 5 x 5 [dgeMatrix] 
#> 
#> Poisson edges: TRUE 
#> Allow self loops: TRUE 
#> 
#> Expected edges: 4995
#> Expected degree: 5
#> Expected density: 0.01
```

Now, conditional on this latent representation, we can sample graphs.
`fastRG` supports several different output types, each of which is
specified by the suffix to `sample_*()` functions. For example, we can
obtain an edgelist in a `tibble` with:

``` r
sample_edgelist(sbm)
#> # A tibble: 4,985 × 2
#>     from    to
#>    <int> <int>
#>  1   111   127
#>  2    86   109
#>  3    43    97
#>  4    61    94
#>  5    22   143
#>  6     4    89
#>  7    30   159
#>  8   119   210
#>  9    41   197
#> 10   145   175
#> # ℹ 4,975 more rows
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
#> IGRAPH 49b52e2 UN-- 1000 5033 -- 
#> + attr: name (v/c)
#> + edges from 49b52e2 (vertex names):
#>  [1] 63 --76  135--215 59 --182 21 --134 180--218 53 --189 138--139 21 --78 
#>  [9] 49 --70  76 --127 6  --139 64 --214 31 --132 56 --93  75 --144 9  --185
#> [17] 33 --150 115--165 163--213 6  --53  47 --179 25 --26  7  --51  10 --55 
#> [25] 120--183 43 --152 25 --34  84 --216 114--191 34 --127 152--164 178--189
#> [33] 106--181 28 --38  41 --89  34 --139 6  --213 24 --153 32 --173 47 --111
#> [41] 157--205 108--133 98 --116 26 --117 18 --194 18 --32  74 --209 18 --128
#> [49] 13 --127 12 --26  1  --133 52 --72  128--213 13 --173 61 --214 33 --142
#> [57] 22 --111 163--191 191--205 5  --108 9  --72  6  --217 113--122 90 --154
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
#> [1]  5.0999835  1.8365365  0.6679806 -0.5241303 -0.8109449
```

Note that eigendecompositions and SVDS (for directed graphs) use
`RSpectra` and do not require explicitly forming large dense population
adjacency matrices; the population decompositions should be efficient in
both time and space for even large graphs.

## Key sampling options

There are several essential tools to modify graph sampling that you
should know about. First there are options that affect the latent factor
sampling:

- `expected_degree`: Set the expected average degree of the graph by
  scaling sampling probabilities. We *strongly, strongly* recommend that
  you always set this option. If you do not, it is easy accidentally
  sample from large and dense graphs.

- `expected_density`: Set the expected density of the graph by scaling
  sampling probabilities. You cannot specify both `expected_degree` and
  `expected_density` at the same time.

In the second stage of graph sampling, the options are:

- `poisson_edges`: Either `TRUE` or `FALSE` depending on whether you
  would like a Bernoulli graph or a Poisson multi-graph. Scaling via
  `expected_degree` assumes a Poisson multi-graph, with some limited
  exceptions.

- `allow_self_edges`: Whether nodes should be allowed to connect to
  themselves. Either `TRUE` or `FALSE`.

## Related work

[`igraph`](https://igraph.org/r/) allows users to sample SBMs (in
$\mathcal O(m + n + k^2)$ time) and random dot product graphs (in
$\mathcal O(n^2 k)$ time).

You can find the original research code associated with `fastRG`
[here](https://github.com/raningtky/sampleRDPG). There is also a Python
translation of the original code in Python
[here](https://github.com/yunjhongwu/matrix-routines/blob/master/fastRG.py).
Both of these implementations are bare bones.

## References

\[1\] Rohe, Karl, Jun Tao, Xintian Han, and Norbert Binkiewicz. 2017. “A
Note on Quickly Sampling a Sparse Matrix with Low Rank Expectation.”
Journal of Machine Learning Research; 19(77):1-13, 2018.
<https://www.jmlr.org/papers/v19/17-128.html>
