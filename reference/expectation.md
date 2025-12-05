# Calculate the expected adjacency matrix

Calculate the expected adjacency matrix

## Usage

``` r
expectation(model, ...)

# S3 method for class 'undirected_factor_model'
expectation(model, ...)

# S3 method for class 'directed_factor_model'
expectation(model, ...)
```

## Arguments

- model:

  A
  [`directed_factor_model()`](https://rohelab.github.io/fastRG/reference/directed_factor_model.md)
  or an
  [`undirected_factor_model()`](https://rohelab.github.io/fastRG/reference/undirected_factor_model.md)
  object.

- ...:

  Unused.

## Value

The expected value of the adjacency matrix, conditional on the latent
factors `X` and `Y` (if the model is directed).
