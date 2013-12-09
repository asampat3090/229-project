% Example file for embedding 10 manifolds from COIL-20 dataset into 2D

% we fix a random seed for repeatability of the results
s = RandStream('mcg16807','Seed',29); RandStream.setGlobalStream(s);

% load coil20 dataset
load coil20

Y0 = double(data);
Y0 = Y0-min(Y0(:));
Y0 = Y0./max(Y0(:));
t0 = objlabel;
M = 10:20; % select a subset of the COIL manifolds
Y = []; t = [];
for i=1:10
  Y = [Y; Y0(t0==M(i),:)];
  t = [t; t0(t0==M(i),:)];
end
N = size(Y,1);

% compute the positive Gaussian affinities with variable width using
% perplexity k=20.
[Wp,beta] = x2p(Y',20);

% compute the negative Gaussian affinities (for EE only) based on the
% square distance matrix
Wn = sqdist(Y);

% Make sure that both matrices are normalized, symmetric and have zeros
% on the diagonal
Wp = (Wp+Wp')/2; Wp(1:N+1:N^2) = 0; Wp = Wp/sum(Wp(:));
Wn = (Wn+Wn')/2; Wn(1:N+1:N^2) = 0; Wn = Wn/sum(Wn(:));

% we want two dimensional embedding
d = 2;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 5 seconds, and the relative tolerance in the embedding is
% no less than 1e-3.
opts.maxit = 100; opts.runtime = 15; opts.tol = 1e-3;

% Dimensionality reduction. Choose between EE, SNE or t-SNE.
l = 100; [X,E,A,T] = ee(Wp,Wn,d,l,opts);
% [X,E,A,T] = ssne(Wp,d,[],opts);
% [X,E,A,T] = tsne(Wp,d,[],opts);

% plot results
figure(1); scatter(X(:,1),X(:,2),24,t);
box on; daspect([1 1 1]);
title(['iter #' num2str(length(E)), ', e=' num2str(E(end)),...
       ', t=' num2str(T(end))]);

