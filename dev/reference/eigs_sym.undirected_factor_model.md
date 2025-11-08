# Compute the eigendecomposition of the expected adjacency matrix of an undirected factor model

Compute the eigendecomposition of the expected adjacency matrix of an
undirected factor model

## Usage

``` r
# S3 method for class 'undirected_factor_model'
eigs_sym(A, k = A$k, which = "LM", sigma = NULL, opts = list(), ...)
```

## Arguments

- A:

  An
  [`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md).

- k:

  Desired rank of decomposition.

- which:

  Selection criterion. See **Details** below.

- sigma:

  Shift parameter. See section **Shift-And-Invert Mode**.

- opts:

  Control parameters related to the computing algorithm. See **Details**
  below.

- ...:

  Unused, included only for consistency with generic signature.

## Details

The `which` argument is a character string that specifies the type of
eigenvalues to be computed. Possible values are:

|      |                                                                                                                                         |
|------|-----------------------------------------------------------------------------------------------------------------------------------------|
| "LM" | The \\k\\ eigenvalues with largest magnitude. Here the magnitude means the Euclidean norm of complex numbers.                           |
| "SM" | The \\k\\ eigenvalues with smallest magnitude.                                                                                          |
| "LR" | The \\k\\ eigenvalues with largest real part.                                                                                           |
| "SR" | The \\k\\ eigenvalues with smallest real part.                                                                                          |
| "LI" | The \\k\\ eigenvalues with largest imaginary part.                                                                                      |
| "SI" | The \\k\\ eigenvalues with smallest imaginary part.                                                                                     |
| "LA" | The \\k\\ largest (algebraic) eigenvalues, considering any negative sign.                                                               |
| "SA" | The \\k\\ smallest (algebraic) eigenvalues, considering any negative sign.                                                              |
| "BE" | Compute \\k\\ eigenvalues, half from each end of the spectrum. When \\k\\ is odd, compute more from the high and then from the low end. |

`eigs()` with matrix types "matrix", "dgeMatrix", "dgCMatrix" and
"dgRMatrix" can use "LM", "SM", "LR", "SR", "LI" and "SI".

[`eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html) with all
supported matrix types, and `eigs()` with symmetric matrix types
("dsyMatrix", "dsCMatrix", and "dsRMatrix") can use "LM", "SM", "LA",
"SA" and "BE".

The `opts` argument is a list that can supply any of the following
parameters:

- `ncv`:

  Number of Lanzcos basis vectors to use. More vectors will result in
  faster convergence, but with greater memory use. For general matrix,
  `ncv` must satisfy \\k+2\le ncv \le n\\, and for symmetric matrix, the
  constraint is \\k \< ncv \le n\\. Default is `min(n, max(2*k+1, 20))`.

- `tol`:

  Precision parameter. Default is 1e-10.

- `maxitr`:

  Maximum number of iterations. Default is 1000.

- `retvec`:

  Whether to compute eigenvectors. If FALSE, only calculate and return
  eigenvalues.

- `initvec`:

  Initial vector of length \\n\\ supplied to the Arnoldi/Lanczos
  iteration. It may speed up the convergence if `initvec` is close to an
  eigenvector of \\A\\.
