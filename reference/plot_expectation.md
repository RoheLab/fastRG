# Plot (expected) adjacency matrices

Plot (expected) adjacency matrices

## Usage

``` r
plot_expectation(model)

plot_dense_matrix(A, ...)

plot_sparse_matrix(A)
```

## Arguments

- model:

  A
  [`directed_factor_model()`](https://rohelab.github.io/fastRG/reference/directed_factor_model.md)
  or an
  [`undirected_factor_model()`](https://rohelab.github.io/fastRG/reference/undirected_factor_model.md)
  object.

- A:

  A [`matrix()`](https://rdrr.io/r/base/matrix.html),
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html) or
  [`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  object.

- ...:

  Unused.

## Value

A
[`ggplot2::ggplot2()`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
plot.

## Examples

``` r
set.seed(27)

model <- dcsbm(n = 10, k = 2, expected_density = 0.2)
#> Generating random degree heterogeneity parameters `theta` from a LogNormal(2, 1) distribution. This distribution may change in the future. Explicitly set `theta` for reproducible results.
#> Generating random mixing matrix `B` with independent Uniform(0, 1) entries. This distribution may change in the future. Explicitly set `B` for reproducible results.

plot_expectation(model)


A <- sample_sparse(model)

plot_sparse_matrix(A)

```
