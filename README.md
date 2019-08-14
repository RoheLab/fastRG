
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastRG

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

`fastRG` quickly samples a broad class of network models known as
generalized random product graphs. In particular, for matrices `X`, `S`
and `Y`, `fastRG` samples a matrix `A` with expectation `X S Y^T` where
individual entries are Poisson distributed. We recommend that you think
of `A` as the adjacency matrix for a graph (or a multi-graph).
Crucially, the sampling is `O(m)`, where `m` is the number of the edges
in graph. Other algorithms are `O(n^2)`, where `n` is the number of
nodes in the network. For additional details, please have a look at the
[paper](https://arxiv.org/abs/1703.02998).

## Installation

`fastRG` is not yet on CRAN. You can install the development version
with:

``` r
# install.package("devtools")
devtools::install_github("RoheLab/fastRG")
```

There’s also a [Python
implementation](https://github.com/yunjhongwu/matrix-routines/blob/master/fastRG.py).

## Example Usage

The easiest way to use `fastRG` is to use wrapper functions that sample
from popular graph models. For example, to sample from an Erdos-Renyi
graph `n = 1,000,000` nodes and expected degree of five, we can use the
`erdos_renyi()` function.

``` r
library(fastRG)
#> Loading required package: Matrix

A <- erdos_renyi(n = 10^6, avg_deg = 5)
```

By default we always get a `Matrix::sparseMatrix()`, but we can also ask
for the graph as an edgelist as well.

``` r
el <- erdos_renyi(n = 1000, avg_deg = 5, return_edge_list = TRUE)
head(el)
#>      from  to
#> [1,]  216 452
#> [2,]  287 571
#> [3,]   13 248
#> [4,]  588 896
#> [5,]  184 208
#> [6,]  418 454
```

This results in a fast way to create `igraph` objects using
`igraph::graph_from_edgelist()`.

``` r
g <- igraph::graph_from_edgelist(el)
g
#> IGRAPH bf08469 D--- 1000 4918 -- 
#> + edges from bf08469:
#>  [1] 216->452 287->571  13->248 588->896 184->208 418->454 440->844
#>  [8]  96->584 174-> 73 760->793  44-> 25 403->763 781->583 226->274
#> [15] 594-> 65  78->871 977->547 212->552 340->831 121->611 917->624
#> [22] 426->676 262->334 867->797 886->119 902->965 181->254 959->519
#> [29]  23->468 361->589  25-> 38 409->684 458->175 108->360 819->684
#> [36] 330->975 714->245 998->784   1->418 671->546 475->898 229->595
#> [43] 785->662 109->463 578->881 971->834 846->898 186->159 218->100
#> [50]  61->330 260->164 436->609 243->745 596-> 15 104->255 195->153
#> [57] 354-> 88 592->223 174->189 768->  8 411->887 720->532 668->124
#> + ... omitted several edges
```
