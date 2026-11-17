# fastRG (development version)

- `sample_edgelist()`, `sample_igraph()` and `sample_tidygraph()` gain a `multigraph` argument. The default, `multigraph = TRUE`, keeps the current behavior: Poisson graphs are returned as multigraphs with repeated edges. With `multigraph = FALSE`, repeated edges are collapsed into a single edge with an integer `weight` attribute (an edgelist column, or an `igraph` edge attribute) that equals the corresponding entry of `sample_sparse()`. In undirected graphs `(i, j)` and `(j, i)` are collapsed together. The `weight` attribute is always present when `multigraph = FALSE`, so Bernoulli graphs (`poisson_edges = FALSE`) get weights of 1 (#30).

- Fixed `sample_sparse()`, `sample_igraph()` and `sample_tidygraph()` erroring on directed rectangular models (`n != d`) when no edges are sampled. Empty graphs now have the same output types as non-empty ones: `sample_edgelist()` returns a tibble with integer `from` and `to` columns instead of a base matrix of doubles, and `sample_sparse()` returns a symmetric `dsCMatrix` for undirected models.

- `overlapping_sbm()` bug fixes:
  - `poisson_edges` and `allow_self_loops` were silently ignored, so every overlapping SBM was a Poisson graph with self-loops. They are now passed on, and Bernoulli overlapping SBMs follow the "noisy-OR" edge probability `1 - prod (1 - B[a, b])^(Z[i, a] Z[j, b])`.
  - Specifying both `k` and `B` is now an error, as in `dcsbm()`.
  - The invertibility error now says that it is the symmetrized and rescaled `B` (i.e. `S`) that must be invertible.
  - `print()` now says "not sorted" for unsorted models, and no longer calls the model degree-corrected.

- `overlapping_sbm()` now puts every node in at least one block. A node that joins no block in the independent Bernoulli draws is assigned to exactly one block, chosen with probability proportional to `pi`. Previously, with the default `pi`, about a third of nodes belonged to no block and had no edges. `pi` must now have at least one positive element. Graphs sampled with a given seed will differ from previous versions.

- New `directed_planted_partition()`, a directed stochastic blockmodel with `k` outgoing and `k` incoming blocks whose mixing matrix has `within_block` (or `a / n`) on the diagonal and `between_block` (or `b / n`) off the diagonal. Outgoing and incoming communities are independent by default; set `shared_blocks = TRUE` to give each node a single community. Supports `block_sizes` / `pi`, `sort_nodes`, `poisson_edges`, `allow_self_loops` and expected degree / density scaling. Like `planted_partition()`, blocks are as balanced as possible by default (#10).

- `planted_partition()` now returns an object with class `undirected_planted_partition` (on top of `undirected_sbm`), so it prints with its own print method, which now shows the within- and between-block edge probabilities. Its validator is now run: `within_block` and `between_block` (or `a / n` and `b / n`) must be scalars between zero and one, as documented, and `B` must have a constant diagonal and off-diagonal. Sampling is unchanged.

- `directed_dcsbm()` gains `shared_blocks` and `shared_theta` arguments, both `FALSE` by default. `shared_blocks = TRUE` samples one block membership vector and uses it for both outgoing and incoming blocks (`z_out` and `z_in` are identical, `B` must be square). `shared_theta = TRUE` uses one degree-heterogeneity vector for both `theta_out` and `theta_in`, and requires `shared_blocks = TRUE`. With shared blocks, `sort_nodes` permutes rows and columns together so that row `i` and column `i` refer to the same node. The defaults keep the co-factor model. The print method now shows whether blocks and degree parameters are shared, and labels outgoing/incoming blocks and degrees correctly (#40).

- The `directed_dcsbm()` documentation now states that in the default co-factor model, `sort_nodes = TRUE` sorts rows and columns independently, so row `i` and column `i` are generally not the same node and self-loops are not meaningful as "node `i` to itself". Use `shared_blocks = TRUE` or `sort_nodes = FALSE` if rows and columns must index the same nodes (#40).

- `directed_dcsbm()` gains `block_sizes_out` and `block_sizes_in` arguments to specify the exact number of nodes in each outgoing and incoming block, as an alternative to sampling block memberships from `pi_out` and `pi_in` (which now default to `NULL`) (#41).

- Breaking change: by default (when neither `pi_out` nor `block_sizes_out` is given, and likewise for the incoming side), `directed_dcsbm()` now makes blocks as balanced as possible, using the same rule as `dcsbm()`, instead of sampling block memberships with equal probabilities. Graphs sampled with a given seed from the `directed_dcsbm()` defaults will differ from previous versions. Passing `pi_out` / `pi_in` explicitly gives the old behavior (#41).

- Fixed `dcsbm()` pairing blocks with the wrong rows and columns of `B` when `sort_nodes = TRUE` (the default). Block memberships were sampled before `B` and `pi` were reordered by block size, and were never relabeled to match. Models with unequal `pi` or `block_sizes`, or with `n` not divisible by `k`, connected the wrong groups of nodes at each rate in `B`. Graphs sampled with a given seed will differ for those models.

- Fixed `dcsbm()` and `directed_dcsbm()` assigning degree parameters to the wrong nodes when `sort_nodes = FALSE`. Theta was written into `X` (and `Y`) in block order rather than node order. With `force_identifiability = FALSE`, nodes got another node's theta; with `force_identifiability = TRUE`, `X` was correct but the returned `theta` was in block order. Sorted models (the default) are unaffected.

- fastRG now requires R >= 4.4.0. irlba 2.4.1, which fastRG uses in tests and vignettes, fails on earlier R versions because `is.atomic(NULL)` changed in R 4.4.0.

- Fixed undersampling of self-loops in undirected models. Since 0.4.0, the diagonal of `A` was sampled at half the intended rate, so `E[A | X, S]` did not equal `X S X'` on the diagonal. `sample_edgelist.undirected_factor_model()` now adds independent Poisson self-loops at the missing rate when `allow_self_loops = TRUE`, so `E[A | X, S] = X S X'` holds exactly. This only affects sampling: `S`, `expected_*()` and population spectra are unchanged, and so are the off-diagonal fixes from #19 and #31. Graphs sampled with a given seed will differ from previous versions (#46).

- Fixed math rendering on the pkgdown site. The home page showed stray backticks because `README.md` used GitHub's `` $`...`$ `` math syntax, which pkgdown can't parse. Also fixed malformed LaTeX (e.g. `\pi_out`, `theta_in`) in the `dcsbm()` and `directed_dcsbm()` documentation (#48).

- Finished the documentation for `overlapping_sbm()`: it now describes the generative model (binary block memberships, additive edge rates), identifiability conditions, related models (OSBM, OCCAM), every argument and the returned fields, and has updated references and examples (#42).

# fastRG 0.4.0

## Major changes

- This release fixes a long-standing inconsistency in how degrees are counted, which resulted in sampling twice as many edges as desired in `undirected_factor_model()`. See below for details (#19). This also caused population singular values for undirected factors models to be off by a factor of 2 (#31).

## Breaking changes

- Specifying both `k` and `B` in `dcsbm()` and `sbm()` now results in an error. You should only specify one of these arguments. 

- It is now possible to specify the precise number of nodes in each block of a `dcsbm()` (and subclasses `sbm()` and `planted_partition()`) via the `block_sizes` argument. This makes it easier to construct blockmodels with exactly repeated eigenvalues. Additionally, the default behavior is now to use this argument and to balance block sizes as evenly as possible. Previously, the default behavior was to sample blocks memberships with equal probability.

## Non-breaking changes

- Added `vignette("consistency")` demonstrating how to check consistency of spectral estimators using `fastRG` for sampling and population spectra computations (#33, #43)

## Details about degree over-sampling bug and the fix

The `fastRG` sampling algorithm, as implemented in `sample_edgelist.matrix()`, is fundamentally a sampler for asymmetric, directed networks with conditional expectation $\mathbb E[A \mid X, S, Y] = X S Y^\top \in \mathbb R^{n_1 \times n_2}$. That is, you can think of the sampler as a very efficient procedure for iterating through $i = 1, ..., n_1$ and $j = 1, ..., n_2$ and sampling from a Poisson with rate $(X S Y^\top)_{ij}$.

However, we would also like to use this same sampler to sample from undirected networks. In an undirected networks, the conditional expectation $\mathbb E[A \mid X, S] = X S X^\top \in \mathbb R^{n \times n}$ is a square matrix with $(X S X^\top)_{ij} = (X S X^\top)_{ji}$. To sample from this matrix, it's typical to sample the upper triangle of $A$ from a Poisson with rate $(X S X^\top)_{ij}$ for all $1 \le i \le j \le n$, and then fill the lower triangle of $A$ to make $A$ symmetric.

Since the `fastRG` algorithm samples $A_{ij}$ for all $i, j$, not just the upper triangle of $A$, we use a trick to sample from undirected networks. First, we force the conditional expectation to be symmetric by symmetrizing $S$. Then, we still sample for all $i, j$. That is, to set $A_{ij}$ we sample once from a Poisson with rate $(X S X^\top)_{ij}$ and once from a Poisson with rate $(X S X^\top)_{ji}$ (these rates are equal by symmetry!). Then we set $A_{ij} = A_{ji}$ to the sum of these Poisson random variables. The issue is that this doubles the expected value of $A_{ij} = A_{ji}$ and so we sample twice as many edges as we should. Up until this release of `fastRG`, we've unfortunately been doing this double sampling in undirected networks (#19).

In this release, we fix this over-sampling. The key is that we divide $S$ by two at sampling time. We do not modify $S$ at all in the `undirected_factor_model()`! You can always use $X S X^\top$ to compute the expected value of $A$. This new change *only affects sampling*. 

That is, instead of passing the $S$ from an `undirected_factor_model()` to the sampler `sample_edgelist.matrix()`, we pass $S / 2$ (see `sample_edgelist.undirected_factor_model()`). This fixes double sampling on the off-diagonal of $A$. The downside is that we're now undersampling by half the diagonal of $A$. I'm assuming that for most use cases this doesn't matter. We could correct for this undersampling of the diagonal of $A$, so please open an issue if self-loops are important to your project.

As a consequence of this change, $A$ and $\mathbb E[A | X, S]$ show now be on the same scale, rather than off by a factor of 2. Importantly, the spectrums should match up now, so you can now use `fastRG` to check how closely you're recovering the spectrum of the your model. See `vignette("consistency")` for a quick demonstration showing consistency of spectral estimates.

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
