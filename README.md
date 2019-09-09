
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastRG

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/RoheLab/fastRG.svg?branch=master)](https://travis-ci.org/RoheLab/fastRG)
<!-- badges: end -->

`fastRG` quickly samples a broad class of network models known as
generalized random product graphs. In particular, for matrices `X`, `S`
and `Y`, `fastRG` samples a matrix `A` with expectation `X S Y^T` where
individual entries are Poisson distributed. We recommend that you think
of `A` as the adjacency matrix for a graph (or a multi-graph).
Crucially, the sampling is `O(m)`, where `m` is the number of the edges
in graph. Other algorithms are `O(n^2)`, where `n` is the number of
nodes in the network. For additional details, see the
[paper](https://arxiv.org/abs/1703.02998).

## Installation

`fastRG` is not yet on CRAN. You can install the development version
with:

``` r
# install.package("devtools")
devtools::install_github("RoheLab/fastRG")
```

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
#> [1,]  143 367
#> [2,]  247 799
#> [3,]  519 224
#> [4,]  334  43
#> [5,]  639 304
#> [6,]  930 702
```

This results in a fast way to create `igraph` objects using
`igraph::graph_from_edgelist()`.

``` r
g <- igraph::graph_from_edgelist(el)
g
#> IGRAPH ef6e76e D--- 1000 5058 -- 
#> + edges from ef6e76e:
#>  [1] 143->367 247->799 519->224 334-> 43 639->304 930->702 673->475
#>  [8]  88->229 723->382 655->631 499-> 40 240->732 710->992  70->463
#> [15] 951->  8 965->755 309->474  67->341 255->654 167->839 610-> 50
#> [22] 163->245 626->635 970->249 349->132 717->870 280->844 186->256
#> [29] 816->320 707->147 400->742 644->585 465->999 126->521 140->986
#> [36] 602->  4 323->311 908->814 183->270  23->330  57->431 365->163
#> [43] 266->489 874->540 397->722 365->502 406-> 30 763->230 772->875
#> [50] 312->  8 467->651 177->  1 839-> 73 818->204  74->650 739->884
#> [57] 520->174 956->621  88->977  95->866 276->924  45->927 969->951
#> + ... omitted several edges
```

## Other languages

The `fastRG` sampler has been implemented in Python
[here](https://github.com/yunjhongwu/matrix-routines/blob/master/fastRG.py).
