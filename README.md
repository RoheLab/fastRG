
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
#> [1,]   57 970
#> [2,]  546 737
#> [3,]  253 779
#> [4,]  798 223
#> [5,]  522 679
#> [6,]  853  72
```

This results in a fast way to create `igraph` objects using
`igraph::graph_from_edgelist()`.

``` r
g <- igraph::graph_from_edgelist(el)
g
#> IGRAPH 6e969ad D--- 1000 5006 -- 
#> + edges from 6e969ad:
#>  [1]  57->970 546->737 253->779 798->223 522->679 853-> 72 574->444
#>  [8] 575->798 946->649 735-> 23  12->729 174->944 121->308  34->968
#> [15] 489->346 231->162 565-> 87 687->320 810->648  64->269 313-> 65
#> [22] 391->701 751->958 488->199 390->375 732->782 870->627 110->396
#> [29] 816->653 421->865 556->907 582->836 938->368 552->602 715->755
#> [36] 144->618 673-> 88 887->825 794->511 710->289 248->243 980->929
#> [43] 868->996 240->925 691->313 454->582 350->345 369->689  67-> 44
#> [50]  17->773 593->367 582->363 983-> 60 472->869 299-> 56 220->713
#> [57] 817->370 770->633  48->479 385->168 338->862  32->750 639->258
#> + ... omitted several edges
```

## Other languages

The `fastRG` sampler has been implemented in Python
[here](https://github.com/yunjhongwu/matrix-routines/blob/master/fastRG.py).
