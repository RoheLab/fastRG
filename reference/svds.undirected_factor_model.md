# Compute the singular value decomposition of the expected adjacency matrix of an undirected factor model

Compute the singular value decomposition of the expected adjacency
matrix of an undirected factor model

## Usage

``` r
# S3 method for class 'undirected_factor_model'
svds(A, k = A$k, nu = k, nv = k, opts = list(), ...)
```

## Arguments

- A:

  An
  [`undirected_factor_model()`](https://rohelab.github.io/fastRG/reference/undirected_factor_model.md).

- k:

  Desired rank of decomposition.

- nu:

  Number of left singular vectors to be computed. This must be between 0
  and `k`.

- nv:

  Number of right singular vectors to be computed. This must be between
  0 and `k`.

- opts:

  Control parameters related to the computing algorithm. See **Details**
  below.

- ...:

  Unused, included only for consistency with generic signature.

## Details

The `opts` argument is a list that can supply any of the following
parameters:

- `ncv`:

  Number of Lanzcos basis vectors to use. More vectors will result in
  faster convergence, but with greater memory use. `ncv` must be satisfy
  \\k \< ncv \le p\\ where `p = min(m, n)`. Default is
  `min(p, max(2*k+1, 20))`.

- `tol`:

  Precision parameter. Default is 1e-10.

- `maxitr`:

  Maximum number of iterations. Default is 1000.

- `center`:

  Either a logical value (`TRUE`/`FALSE`), or a numeric vector of length
  \\n\\. If a vector \\c\\ is supplied, then SVD is computed on the
  matrix \\A - 1c'\\, in an implicit way without actually forming this
  matrix. `center = TRUE` has the same effect as `center = colMeans(A)`.
  Default is `FALSE`.

- `scale`:

  Either a logical value (`TRUE`/`FALSE`), or a numeric vector of length
  \\n\\. If a vector \\s\\ is supplied, then SVD is computed on the
  matrix \\(A - 1c')S\\, where \\c\\ is the centering vector and \\S =
  diag(1/s)\\. If `scale = TRUE`, then the vector \\s\\ is computed as
  the column norm of \\A - 1c'\\. Default is `FALSE`.
