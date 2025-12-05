# fastRG

`fastRG` quickly samples a broad class of network models known as
generalized random dot product graphs (GRDPGs). In particular, for
matrices \$\`X\`\$, \$\`S\`\$ and \$\`Y\`\$, `fastRG` samples a matrix
\$\`A\`\$ with expectation \$\`X S Y^T\`\$ where the entries are
independently Poisson distributed conditional on \$\`X\`\$ and
\$\`Y\`\$. This is primarily useful when \$\`A\`\$ is the adjacency
matrix of a graph. Crucially, the sampling is \$\`\mathcal O(m)\`\$,
where \$\`m\`\$ is the number of the edges in graph, as opposed to the
naive sampling approach, which is \$\`\mathcal O(n^2)\`\$, where
\$\`n\`\$ is the number of nodes in the network. For additional details,
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
graphs. First, we sample the latent factors \$\`X\`\$ and \$\`Y\`\$.
Then we sample \$\`A\`\$ conditional on those latent factors. `fastRG`
mimics this two-stage sample structure. For example, to sample from a
stochastic blockmodel, we first create the latent factors.

``` r
library(fastRG)
#> Loading required package: Matrix

set.seed(27)

sbm <- sbm(n = 1000, k = 5, expected_density = 0.01)
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.
```

You can specify the latent factors and the mixing matrix \$\`B\`\$
yourself, but there are also defaults to enable fast prototyping. Here
\$\`B\`\$ was randomly generated with `Uniform[0, 1]` entries and nodes
were assigned randomly to communities with equal probability of falling
in all communities. Printing the result object gives us some additional
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
#> # A tibble: 2,484 × 2
#>     from    to
#>    <int> <int>
#>  1     4   155
#>  2    46   141
#>  3    42    56
#>  4    42    55
#>  5    72   167
#>  6    32    68
#>  7    67    75
#>  8    10   164
#>  9    30   154
#> 10    74   182
#> # ℹ 2,474 more rows
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
#> IGRAPH 2e1d765 UN-- 1000 2447 -- 
#> + attr: name (v/c)
#> + edges from 2e1d765 (vertex names):
#>  [1] 125--139 30 --36  36 --57  76 --108 46 --128 44 --179 49 --199 9  --69 
#>  [9] 17 --154 100--154 58 --138 22 --182 27 --92  109--143 44 --195 96 --153
#> [17] 7  --68  121--159 17 --42  132--171 53 --145 15 --33  68 --78  58 --99 
#> [25] 158--169 34 --159 128--194 23 --74  6  --126 33 --139 33 --128 80 --107
#> [33] 8  --55  45 --156 120--133 8  --88  120--138 15 --26  123--173 26 --68 
#> [41] 145--148 77 --123 1  --110 20 --41  90 --184 72 --191 37 --90  36 --192
#> [49] 101--119 116--131 159--188 37 --58  50 --170 6  --40  132--154 157--194
#> [57] 130--136 14 --143 89 --195 143--173 72 --81  30 --184 159--176 34 --126
#> + ... omitted several edges
```

Note that every time we call `sample_*()` we draw a new sample.

``` r
A <- sample_sparse(sbm)
B <- sample_sparse(sbm)

all(A == B) # random realizations from the SBM don't match!
#> [1] FALSE
```

## Efficient spectral decompositions

If you would like to obtain the singular value decomposition of the
population adjacency matrix conditional on latent factors, that is
straightforward:

``` r
s <- eigs_sym(sbm)
s$values
#> [1]  5.0838486  1.8176036  0.6987030 -0.5157282 -0.8208442
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
\$\`\mathcal O(m + n + k^2)\`\$ time) and random dot product graphs (in
\$\`\mathcal O(n^2 k)\`\$ time).

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
