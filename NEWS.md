# fastRG (development version)

## Major changes

- This release fixes a long-standing inconsistency in how degrees are counted, which resulted in sampling twice as many edges as desired in `undirected_factor_model()`. See below for details (#19). This also caused population singular values for undirected factors models to be off by a factor of 2 (#31).

## Breaking changes

- Specifying both `k` and `B` in `dcsbm()` and `sbm()` now results in an error. You should only specify one of these arguments. 

- It is now possible to specify the precise number of nodes in each block of a `dcsbm()` (and subclasses `sbm()` and `planted_partition()`) via the `block_sizes` argument. This makes it easier to construct blockmodels with exactly repeated eigenvalues. Additionally, the default behavior is now to use this argument and to balance block sizes as evenly as possible. Previously, the default behavior was to sample blocks memberships with equal probability.

## Non-breaking changes

- Added `vignette("consistency")` demonstratingn how to check consistency of spectral estimators using `fastRG` for sampling and population spectra computations (#33, #43)

## Details about degree over-sampling bug and the fix

The `fastRG` sampling algorithm, as implemented in `sample_edgelist.matrix()`, is fundamentally a sampler for asymmetric, directed networks with conditional expectation $\mathbb E[A \mid X, S, Y] = X S Y^top \in \mathbb R^{n_1 \times n_2}$. That is, you can think of the sampler as a very efficient procedure for iterating through $i = 1, ..., n_1$ and $j = 1, ..., n_2$ and sampling from a Poisson with rate $(X S Y^\top)_{ij}$.

However, we would also like to use this same sampler to sample from undirected networks. In an undirected networks, we have the conditional expectation $\mathbb E[A \mid X, S, Y] = Z B Y^top \in \mathbb R^{n \times n}$ is a square matrix with $(X S Y^\top)_{ij} = (X S Y^\top)_{ji}$. To sample from this matrix, you want to sample the upper triangle of $A$ from a Poisson with rate $(X S Y^\top)_{ij}$ for all $1 \le i \le j \le n$, and then symmetrize $A$.

Since the `fastRG` algorithm samples $A_{ij}$ for all $i, j$, not just the upper triangle of $A$, we use a trick to sample from undirected networks. First, we force the conditional expectation to the symmetric by symmetrizing $S$. Then, we still sample for all $i, j$. To set $A_{ij}$ we sample once from a Poisson with rate $(X S Y^\top)_{ij}$ and once from a Poisson with rate $(X S Y^\top)_{ji}$ (these rates are equal by symmetry!). Then we set $A_{ij} = A_{ji}$ to the sum of these sample Poisson random variables. The issue is that this doubles the expected value of $A_{ij} = A_{ji}$ and so we sample twice as many edges as we should. Up until this release of `fastRG`, we've unfortunately been doing this double sampling (see #19).

In this release, we fix this over-sampling. The key is that we divide $S$ by two at sampling time. We do not modify $S$ at all in the `undirected_factor_model()`! You can always use $X S X^\top$ to compute the expected value of $A$. This new change *only affects sampling*. 

That is, instead of passing the $S$ from an `undirected_factor_model()` to the sampler `sample_edgelist.matrix()`, we pass $S / 2$ (see `sample_edgelist.undirected_factor_model()`). This fixes double sampling on the off-diagonal of $A$. The downside is that we're now undersampling by half the diagonal of $A$. I'm assuming that for most use cases this doesn't matter. We could correct for this undersampling of the diagonal of $A$, so please open an issue if self-loops are important to your project.

As a consequence of this change, $A$ and $\mathbb E[A | X, S, Y]$ show now be on the same scale, rather than off by a factor of 2. Importantly, the spectrums should match up now, so you can now use `fastRG` to check how closely you're recovering the spectrum of the your model. See `vignette("consistency")` for a quick demonstration show consistency of spectral estimates.

# fastRG 0.3.3

- Improve cross-linking to documentation of other packages for CRAN

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
