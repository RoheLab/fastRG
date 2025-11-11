# Package index

## Sampling

- [`sample_edgelist()`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.md)
  : Sample a random edgelist from a random dot product graph
- [`sample_edgelist(`*`<matrix>`*`)`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.matrix.md)
  [`sample_edgelist(`*`<Matrix>`*`)`](https://rohelab.github.io/fastRG/dev/reference/sample_edgelist.matrix.md)
  : Low level interface to sample RPDG edgelists
- [`sample_igraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_igraph.md)
  : Sample a random dot product graph as an igraph graph
- [`sample_sparse()`](https://rohelab.github.io/fastRG/dev/reference/sample_sparse.md)
  : Sample a random dot product graph as a sparse Matrix
- [`sample_tidygraph()`](https://rohelab.github.io/fastRG/dev/reference/sample_tidygraph.md)
  : Sample a random dot product graph as a tidygraph graph

## Graph properties conditional on latent representations

- [`eigs_sym(`*`<undirected_factor_model>`*`)`](https://rohelab.github.io/fastRG/dev/reference/eigs_sym.undirected_factor_model.md)
  : Compute the eigendecomposition of the expected adjacency matrix of
  an undirected factor model
- [`svds(`*`<directed_factor_model>`*`)`](https://rohelab.github.io/fastRG/dev/reference/svds.directed_factor_model.md)
  : Compute the singular value decomposition of the expected adjacency
  matrix of a directed factor model
- [`svds(`*`<undirected_factor_model>`*`)`](https://rohelab.github.io/fastRG/dev/reference/svds.undirected_factor_model.md)
  : Compute the singular value decomposition of the expected adjacency
  matrix of an undirected factor model
- [`expected_edges()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md)
  [`expected_degrees()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md)
  [`expected_degree()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md)
  [`expected_in_degree()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md)
  [`expected_out_degree()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md)
  [`expected_density()`](https://rohelab.github.io/fastRG/dev/reference/expected_edges.md)
  : Calculate the expected edges in Poisson RDPG graph
- [`expectation()`](https://rohelab.github.io/fastRG/dev/reference/expectation.md)
  : Calculate the expected adjacency matrix
- [`plot_expectation()`](https://rohelab.github.io/fastRG/dev/reference/plot_expectation.md)
  [`plot_dense_matrix()`](https://rohelab.github.io/fastRG/dev/reference/plot_expectation.md)
  [`plot_sparse_matrix()`](https://rohelab.github.io/fastRG/dev/reference/plot_expectation.md)
  : Plot (expected) adjacency matrices

## Undirected graphs

- [`erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/erdos_renyi.md)
  : Create an undirected erdos renyi object
- [`chung_lu()`](https://rohelab.github.io/fastRG/dev/reference/chung_lu.md)
  : Create an undirected Chung-Lu object
- [`planted_partition()`](https://rohelab.github.io/fastRG/dev/reference/planted_partition.md)
  : Create an undirected planted partition object
- [`sbm()`](https://rohelab.github.io/fastRG/dev/reference/sbm.md) :
  Create an undirected stochastic blockmodel object
- [`dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/dcsbm.md) :
  Create an undirected degree corrected stochastic blockmodel object
- [`overlapping_sbm()`](https://rohelab.github.io/fastRG/dev/reference/overlapping_sbm.md)
  : Create an undirected overlapping degree corrected stochastic
  blockmodel object
- [`undirected_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/undirected_factor_model.md)
  : Create an undirected factor model graph
- [`mmsbm()`](https://rohelab.github.io/fastRG/dev/reference/mmsbm.md) :
  Create an undirected degree-corrected mixed membership stochastic
  blockmodel object

## Directed graphs

- [`directed_erdos_renyi()`](https://rohelab.github.io/fastRG/dev/reference/directed_erdos_renyi.md)
  : Create an directed erdos renyi object
- [`directed_dcsbm()`](https://rohelab.github.io/fastRG/dev/reference/directed_dcsbm.md)
  : Create a directed degree corrected stochastic blockmodel object
- [`directed_factor_model()`](https://rohelab.github.io/fastRG/dev/reference/directed_factor_model.md)
  : Create a directed factor model graph
