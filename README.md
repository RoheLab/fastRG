
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

We *strongly* advise that you always set `avg_deg`, as it is easy to
request very large and dense graphs without this scaling.

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
#> [1,]  915 640
#> [2,]  606 727
#> [3,]  957 530
#> [4,]  420 609
#> [5,]  952 138
#> [6,]  455 566
```

This results in a fast way to create `igraph` objects using
`igraph::graph_from_edgelist()`.

``` r
g <- igraph::graph_from_edgelist(el)
g
#> IGRAPH 2123d20 D--- 1000 5116 -- 
#> + edges from 2123d20:
#>  [1] 915->640 606->727 957->530 420->609 952->138 455->566 228->908
#>  [8] 781->465 393->126 905->358 331-> 96 783->541 393->678 110->493
#> [15] 310->293   9-> 31 679->542 510->836 916->718 923->255 566->820
#> [22] 726-> 53  34->368 676->998 321->476 963->281 183->408 988->258
#> [29] 461-> 49 316->630 167->698  45->244 840->697 840->780 523->479
#> [36] 176->859  42->423 997->162  57-> 53 118->969 490->784 193->266
#> [43] 728->283 926->411 140->348 557->254 222->664 431->910 612-> 54
#> [50] 808->534 774->213 974->497 798->878 441->766  81->980 896->846
#> [57] 510->501 864->863 708->133 352->112 947->794 978->983 436->189
#> + ... omitted several edges
```
