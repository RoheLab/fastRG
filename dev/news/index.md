# Changelog

## fastRG (development version)

- Added option to specify precise number of nodes in each block of a
  [`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md)
  or [`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)
  via the `block_sizes` argument. This makes it easier to construct
  blockmodels with exactly repeated eigenvalues.
- The default behavior of
  [`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md),
  [`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md) and
  [`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md)
  has changed: when `block_sizes` or `pi` is unspecified, the new
  default is to balance block sizes as evenly as possible. Previously,
  `pi` was set to a constant vector, balancing block sizes in
  expectation only.
- Specifying both `k` and `B` in
  [`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md)
  and [`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md)
  now results in an error; only specify one of these arguments.

## fastRG 0.3.3

CRAN release: 2025-07-24

- Improve cross-linking to documentation of other packages for CRAN

## fastRG 0.3.2

CRAN release: 2023-08-21

- Added documentation about block sorting in blockmodels when
  `sort_nodes = TRUE`
  ([\#35](https://github.com/RoheLab/fastRG/issues/35)). Blocks are now
  only sorted when `sort_nodes = TRUE`, although they were previously
  always sorted. In directed stochastic blocks, flipped incoming and
  outgoing blocks, such that `X` now contains info about outgoing blocks
  and `Y` now contains info about incoming blocks, as you would expected
  if `A[i, j]` encodes an edge from node `i` to node `j`
- Fixed bug where isolated nodes were sometimes dropped from igraph and
  tidygraph objects
  ([\#35](https://github.com/RoheLab/fastRG/issues/35))
- Added
  [`plot_expectation()`](https://rohelab.github.io/fastRG/dev/reference/plot_expectation.md),
  [`plot_sparse_matrix()`](https://rohelab.github.io/fastRG/dev/reference/plot_expectation.md)
  and
  [`expectation()`](https://rohelab.github.io/fastRG/dev/reference/expectation.md)
  utilities ([\#34](https://github.com/RoheLab/fastRG/issues/34))
- Fixed incorrect computation in
  [`expected_degrees()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md)
  ([\#34](https://github.com/RoheLab/fastRG/issues/34))

## fastRG 0.3.1

CRAN release: 2022-06-30

### Breaking changes

- Users must now pass `poisson_edges` and `allow_self_loops` arguments
  to model object constructors
  (i.e. [`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md))
  rather than `sample_*()` methods. Additionally, when
  `poisson_edges = FALSE`, the mixing matrix `S` is taken (after
  degree-scaling and possible symmetrization for undirected models) to
  represent desired inter-factor connection probabilities, and thus
  should be between zero and one. This Bernoulli-parameterized `S` is
  then transformed into the equivalent (or approximately equivalent)
  Poisson `S`. See Section 2.3 of Rohe et al. (2017) for additional
  details about this conversion and approximation of Bernoulli graphs by
  Poisson graphs ([\#29](https://github.com/RoheLab/fastRG/issues/29)).

### Other news

- Add overlapping stochastic blockmodel
  ([\#7](https://github.com/RoheLab/fastRG/issues/7),
  [\#25](https://github.com/RoheLab/fastRG/issues/25))
- Add directed degree-corrected stochastic blockmodels
  ([\#18](https://github.com/RoheLab/fastRG/issues/18))
- Allow rank 1 undirected stochastic block models
- Fix bug where isolated nodes where dropped from sampled tidygraphs
  ([\#23](https://github.com/RoheLab/fastRG/issues/23))
- Allow users to force model identification in DC-SBMs with
  `force_identifiability = TRUE`, and in overlapping SBMs with
  `force_pure = TRUE`, which are now the default.
- Improve population expected degree/density computations
  ([\#19](https://github.com/RoheLab/fastRG/issues/19))
- Let user know when `theta_out` is automatically generated for directed
  DC-SBMs ([\#22](https://github.com/RoheLab/fastRG/issues/22))
- Fixed an obscure but pesky issue sampling from models with empty
  blocks ([\#13](https://github.com/RoheLab/fastRG/issues/13))
- Documented [`svds()`](https://rdrr.io/pkg/RSpectra/man/svds.html) and
  [`eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html) methods,
  which allow users to take spectral decompositions of expected
  adjacency matrices conditional `X`, `S` and `Y`.

## fastRG 0.3.0

CRAN release: 2021-02-26

- Released to CRAN

## fastRG 0.2.0.9000

- Added a `NEWS.md` file to track changes to the package.
