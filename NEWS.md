# fastRG (development version)

* Add overlapping stochastic blockmodel (#7, #25)
* Add directed degree-corrected stochastic blockmodels (#18)
* Allow rank 1 undirected stochastic block models
* Fix bug where isolated nodes where dropped from sampled tidygraphs (#23)
* Allow users to force model identification in DC-SBMs with `force_identifiability = TRUE`, and in overlapping SBMs with `force_pure = TRUE`, which are now the default.
* Improve population expected degree/density computations (#19)
* Let user know when `theta_out` is automatically generated for directed DC-SBMs (#22)
* Document `svds()` and `eigs_sym()` methods, which allow users to take spectral decompositions of expected adjacency matrices conditional `X`, `S` and `Y`.

# fastRG 0.3.0

* Released to CRAN

# fastRG 0.2.0.9000

* Added a `NEWS.md` file to track changes to the package.
