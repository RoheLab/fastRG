
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
from population graph models. For example, to sample from an Erdos-Renyi
graph `n = 1,000,000` nodes and expected degree 5, we can use the `er()`
function.

``` r
library(fastRG)
#> Loading required package: Matrix

A <- erdos_renyi(n = 10^6, avgDeg = 5)
```

By default we always get a `Matrix::sparseMatrix()`, but we can also ask
for the graph as an edgelist as well. This results in a fast way to
create `igraph` objects using `igraph::graph_from_edgelist()`.

``` r
el <- erdos_renyi(n = 1000, avgDeg = 5, returnEdgeList = TRUE)
head(el)
#>       eo  ei
#> [1,] 596 926
#> [2,] 108 449
#> [3,]  58 313
#> [4,] 444 852
#> [5,] 296 548
#> [6,] 592 273
```

``` r
g <- igraph::graph_from_edgelist(el)
g
#> IGRAPH 717e098 D--- 1000 4986 -- 
#> + edges from 717e098:
#>  [1] 596->926 108->449  58->313 444->852 296->548 592->273 920->904
#>  [8]  80->621 472->697  79->186 680->180 763->905 184->865 776->863
#> [15] 842->903 764->268  81->420 725->433 604->521 377->902 524->921
#> [22]  23->480 152->869  51->762 299->128 372->572 188->546 925->777
#> [29] 278->630 735->213 841->186  68->147 426->263 529->654 802->450
#> [36] 302->874 629->725 267->846 804->794  57->397 658-> 98 983->533
#> [43] 257-> 11 648->956 313->592 392-> 26 834->729 403-> 21 368->577
#> [50] 924->286 676->114 141->164 193->488 379->570 670->761 250-> 41
#> [57]  77->389  44->289  33->788 632-> 26 570->815 451->862 637->625
#> + ... omitted several edges
```
