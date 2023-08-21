# fastRG 0.3.2

- Added documentation about block sorting in blockmodels when `sort_nodes = TRUE` (#35). Blocks are now only sorted when `sort_nodes = TRUE`, although they were previously always sorted. In directed stochastic blocks, flipped incoming and outgoing blocks, such that `X` now contains info about outgoing blocks and `Y` now contains info about incoming blocks, as you would expected if `A[i, j]` encodes an edge from node `i` to node `j`
- Fixed bug where isolated nodes were sometimes dropped from igraph and tidygraph objects (#35)
- Added `plot_expectation()`, `plot_sparse_matrix()` and `expectation()` utilities (#34)
- Fixed incorrect computation in `expected_degrees()` (#34)

# fastRG 0.3.1

## Breaking changes

- Users must now pass `poisson_edges` and `allow_self_loops` arguments to model object constructors (i.e. `sbm()`) rather than `sample_*()` methods. Additionally, when `poisson_edges = FALSE`, the mixing matrix `S` is taken (after degree-scaling and possible symmetrization for undirected models) to represent desired inter-factor connection probabilities, and thus should be between zero and one. This Bernoulli-parameterized `S` is then transformed into the equivalent (or approximately equivalent) Poisson `S`. See Section 2.3 of Rohe et al. (2017) for additional details about this conversion and approximation of Bernoulli graphs by Poisson graphs (#29).

## Other news

* Add overlapping stochastic blockmodel (#7, #25)
* Add directed degree-corrected stochastic blockmodels (#18)
* Allow rank 1 undirected stochastic block models
* Fix bug where isolated nodes where dropped from sampled tidygraphs (#23)
* Allow users to force model identification in DC-SBMs with `force_identifiability = TRUE`, and in overlapping SBMs with `force_pure = TRUE`, which are now the default.
* Improve population expected degree/density computations (#19)
* Let user know when `theta_out` is automatically generated for directed DC-SBMs (#22)
* Fixed an obscure but pesky issue sampling from models with empty blocks (#13)
* Documented `svds()` and `eigs_sym()` methods, which allow users to take spectral decompositions of expected adjacency matrices conditional `X`, `S` and `Y`.

# fastRG 0.3.0

* Released to CRAN

# fastRG 0.2.0.9000

* Added a `NEWS.md` file to track changes to the package.
