# Create an undirected stochastic blockmodel object

To specify a stochastic blockmodel, you must specify the number of nodes
(via `n`), the mixing matrix (via `k` or `B`), and the relative block
probabilites (optional, via `pi`). We provide defaults for most of these
options to enable rapid exploration, or you can invest the effort for
more control over the model parameters. We **strongly recommend**
setting the `expected_degree` or `expected_density` argument to avoid
large memory allocations associated with sampling large, dense graphs.

## Usage

``` r
sbm(
  n,
  k = NULL,
  B = NULL,
  ...,
  block_sizes = NULL,
  pi = NULL,
  sort_nodes = TRUE,
  poisson_edges = TRUE,
  allow_self_loops = TRUE
)
```

## Arguments

- n:

  The number of nodes in the network. Must be a positive integer. This
  argument is required.

- k:

  (mixing matrix) The number of blocks in the blockmodel. Use when you
  don't want to specify the mixing-matrix by hand. When `k` is
  specified, the elements of `B` are drawn randomly from a
  `Uniform(0, 1)` distribution. This is subject to change, and may not
  be reproducible. `k` defaults to `NULL`. You must specify either `k`
  or `B`, but not both.

- B:

  (mixing matrix) A `k` by `k` matrix of block connection probabilities.
  The probability that a node in block `i` connects to a node in
  community `j` is `Poisson(B[i, j])`. Must be a square matrix. `matrix`
  and `Matrix` objects are both acceptable. If `B` is not symmetric, it
  will be symmetrized via the update `B := B + t(B) / 2`. Defaults to
  `NULL`. You must specify either `k` or `B`, but not both.

- ...:

  Arguments passed on to
  [`undirected_factor_model`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)

  `expected_degree`

  :   If specified, the desired expected degree of the graph. Specifying
      `expected_degree` simply rescales `S` to achieve this. Defaults to
      `NULL`. Do not specify both `expected_degree` and
      `expected_density` at the same time.

  `expected_density`

  :   If specified, the desired expected density of the graph.
      Specifying `expected_density` simply rescales `S` to achieve this.
      Defaults to `NULL`. Do not specify both `expected_degree` and
      `expected_density` at the same time.

- block_sizes:

  (block sizes) Number of nodes in each block, as a vector of integers.
  Must match the dimensions of `B`, or `k` and must sum to `n`. Defaults
  to `NULL`, in which case blocks are made to be as balanced as
  possible. You can specify either `pi` or `block_sizes`, but not both.

- pi:

  (block sizes) Relative block probabilities. Must be positive, but do
  not need to sum to one, as they will be normalized internally. Must
  match the dimensions of `B` or `k`. Defaults to `NULL`, in which case
  the `block_sizes` argument will take precedence. Note that you can
  specify either `pi` or `block_sizes`, but should not specify both.

- sort_nodes:

  Logical indicating whether or not to sort the nodes so that they are
  grouped by block and by `theta`. Useful for plotting. Defaults to
  `TRUE`. When `TRUE`, nodes are first sorted by block membership, and
  then by degree-correction parameters within each block. Additionally,
  `pi` is sorted in increasing order, and the columns of the `B` matrix
  are permuted to match the new order of `pi`.

- poisson_edges:

  Logical indicating whether or not multiple edges are allowed to form
  between a pair of nodes. Defaults to `TRUE`. When `FALSE`, sampling
  proceeds as usual, and duplicate edges are removed afterwards.
  Further, when `FALSE`, we assume that `S` specifies a desired
  between-factor connection probability, and back-transform this `S` to
  the appropriate Poisson intensity parameter to approximate Bernoulli
  factor connection probabilities. See Section 2.3 of Rohe et al. (2017)
  for some additional details.

- allow_self_loops:

  Logical indicating whether or not nodes should be allowed to form
  edges with themselves. Defaults to `TRUE`. When `FALSE`, sampling
  proceeds allowing self-loops, and these are then removed after the
  fact.

## Value

An `undirected_sbm` S3 object, which is a subclass of the
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md)
object.

## Details

A stochastic block is equivalent to a degree-corrected stochastic
blockmodel where the degree heterogeneity parameters have all been set
equal to 1.

## See also

Other stochastic block models:
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md)

Other undirected graphs:
[`chung_lu()`](https://rohelab.github.io/fastRG/dev/reference/chung_lu.md),
[`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
[`erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/erdos_renyi.md),
[`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md),
[`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md),
[`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md)

## Examples

``` r
set.seed(27)

lazy_sbm <- sbm(n = 100, k = 5, expected_density = 0.01)
lazy_sbm
#> $fun
#> function (n, pref.matrix, block.sizes, directed = FALSE, loops = FALSE) 
#> {
#>     n <- as.numeric(n)
#>     pref.matrix[] <- as.numeric(pref.matrix)
#>     block.sizes <- as.numeric(block.sizes)
#>     directed <- as.logical(directed)
#>     loops <- as.logical(loops)
#>     on.exit(.Call(R_igraph_finalizer))
#>     res <- .Call(R_igraph_sbm_game, n, pref.matrix, block.sizes, 
#>         directed, loops)
#>     if (igraph_opt("add.params")) {
#>         res$name <- "Stochastic block model"
#>         res$loops <- loops
#>     }
#>     res
#> }
#> <bytecode: 0x556822c1c320>
#> <environment: namespace:igraph>
#> 
#> $args
#> <list_of<quosure>>
#> 
#> $n
#> <quosure>
#> expr: ^100
#> env:  empty
#> 
#> $k
#> <quosure>
#> expr: ^5
#> env:  empty
#> 
#> $expected_density
#> <quosure>
#> expr: ^0.01
#> env:  empty
#> 
#> 
#> $lazy
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "igraph_constructor_spec"

# by default we get a multigraph (i.e. multiple edges are
# allowed between the same two nodes). using bernoulli edges
# will with an adjacency matrix with only zeroes and ones

bernoulli_sbm <- sbm(
  n = 500,
  k = 30,
  poisson_edges = FALSE,
  expected_degree = 8
)

bernoulli_sbm
#> $fun
#> function (n, pref.matrix, block.sizes, directed = FALSE, loops = FALSE) 
#> {
#>     n <- as.numeric(n)
#>     pref.matrix[] <- as.numeric(pref.matrix)
#>     block.sizes <- as.numeric(block.sizes)
#>     directed <- as.logical(directed)
#>     loops <- as.logical(loops)
#>     on.exit(.Call(R_igraph_finalizer))
#>     res <- .Call(R_igraph_sbm_game, n, pref.matrix, block.sizes, 
#>         directed, loops)
#>     if (igraph_opt("add.params")) {
#>         res$name <- "Stochastic block model"
#>         res$loops <- loops
#>     }
#>     res
#> }
#> <bytecode: 0x556822c1c320>
#> <environment: namespace:igraph>
#> 
#> $args
#> <list_of<quosure>>
#> 
#> $n
#> <quosure>
#> expr: ^500
#> env:  empty
#> 
#> $k
#> <quosure>
#> expr: ^30
#> env:  empty
#> 
#> $poisson_edges
#> <quosure>
#> expr: ^FALSE
#> env:  empty
#> 
#> $expected_degree
#> <quosure>
#> expr: ^8
#> env:  empty
#> 
#> 
#> $lazy
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "igraph_constructor_spec"

edgelist <- sample_edgelist(bernoulli_sbm)
#> Error in UseMethod("sample_edgelist"): no applicable method for 'sample_edgelist' applied to an object of class "igraph_constructor_spec"
edgelist
#> Error: object 'edgelist' not found

A <- sample_sparse(bernoulli_sbm)
#> Error in UseMethod("sample_sparse"): no applicable method for 'sample_sparse' applied to an object of class "igraph_constructor_spec"

# only zeroes and ones!
sign(A)
#> Error: object 'A' not found


# sbm with repeated eigenvalues

# block sizes equal by default, needed to prevent variation in spectrum
# from variation in block sizes. also need B to have a single repeated
# eigenvalue

repeated_eigen <- sbm(
  n = 100,
  B = diag(rep(0.8, 5)),
  expected_degree = 10
)

# exactly repeated eigenvalues in the population
e <- eigs_sym(repeated_eigen)
#> Error in UseMethod("eigs_sym"): no applicable method for 'eigs_sym' applied to an object of class "igraph_constructor_spec"
e$values
#> Error: object 'e' not found
```
