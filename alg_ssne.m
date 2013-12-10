% Computes tsne using Max Vladymyrov's fast tsne algorithm
% assumes data is an (n x m) matrix with n data points in m-dimensional
% space
function [X, E, A, T] = alg_ssne(data, dim, opts)

display('SSNE: pre-process data');

% compute the positive Gaussian affinities with k=20 nearest neighbors and
% sigma=5
Wp = gaussaff(data,{'K',20},5);

% Make sure that both matrices are normalized, symmetric and have zeros
% on the diagonal
N = size(data, 1);
Wp = (Wp+Wp')/2; Wp(1:N+1:N^2) = 0; Wp = Wp/sum(Wp(:));

% Dimensionality reduction.
display('SSNE: Dimensionality Reduction');
[X,E,A,T] = ssne(Wp,dim,[],opts);
display('SSNE: Finished');