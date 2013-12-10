% Computes tsne using Max Vladymyrov's fast tsne algorithm
% assumes data is an (n x m) matrix with n data points in m-dimensional
% space
function [X, E, A, T] = alg_ee(data, dim, opts)

display('EE: pre-process data');

% compute the positive Gaussian affinities with k=20 nearest neighbors and
% sigma=5
Wp = gaussaff(data,{'K',20},5);

% compute the negative Gaussian affinities (for EE only) based on the
% square distance matrix
Wn = sqdist(data);

% Make sure that both matrices are normalized, symmetric and have zeros
% on the diagonal
N = size(data, 1);
Wp = (Wp+Wp')/2; Wp(1:N+1:N^2) = 0; Wp = Wp/sum(Wp(:));
Wn = (Wn+Wn')/2; Wn(1:N+1:N^2) = 0; Wn = Wn/sum(Wn(:));

% Dimensionality reduction.
display('EE: Dimensionality Reduction');
l = 1e5; 
[X,E,A,T] = ee(Wp,Wn,dim,l,opts);
display('EE: Finished');