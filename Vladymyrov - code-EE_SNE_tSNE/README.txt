(C) 2012 by Miguel A. Carreira-Perpinan and Max Vladymyrov
   Electrical Engineering and Computer Science
   University of California, Merced
   http://eecs.ucmerced.edu

The functions listed below perform a fast optimization for the following
nonlinear embedding algorithms:
- EE: elastic embedding
- s-SNE: symmetric stochastic neighbor embedding
- t-SNE: t-distributed stochastic neighbor embedding

The fast optimization is based on the spectral direction described in this
paper, with backtracking line search:

  Max Vladymyrov and Miguel Á. Carreira-Perpiñán: "Partial-Hessian
  strategies for fast learning of nonlinear embeddings". ICML 2012.

This code does not implement other variations of the spectral direction
described in the paper (e.g. SD-).

See each function for detailed usage instructions. The functions require
as input an affinity matrix Wp of NxN weights obtained from the N
high-dimensional data points (not the points themselves). The examples
show how to construct the matrix Wp, and how to run EE, s-SNE or t-SNE.

List of functions:
- ee.m, ssne.m tsne.m: training functions for EE, s-SNE and t-SNE.
- demo*: examples using spiral, Swiss roll and COIL-20 datasets [1].
  The latter requires the file coil20.mat available here:
  http://faculty.ucmerced.edu/mcarreira-perpinan/software/coil20.mat

The following are used internally by other functions:
- sqdist: matrix of Euclidean squared distances between two point sets.
- gaussaff: computes the (sparse) affinity matrix from the given dataset.
- x2p.m: computes an adaptive bandwidth of the Gaussian kernel for a given
  perplexity value. (Code originally from Sam Roweis.)

[1] S. A. Nene, S. K. Nayar and H. Murase: "Columbia Object Image Library
    (COIL-20)", Technical Report CUCS-005-96, Feb. 1996.
    http://www.cs.columbia.edu/CAVE/software/softlib/coil-20.php

