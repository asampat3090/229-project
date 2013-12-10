% Example file for embedding 3D Swiss roll dataset into 2D

% we fix a random seed for repeatability of the results
s = RandStream('mcg16807','Seed',29); 
RandStream.setGlobalStream(s);
% RandStream.setDefaultStream(s); % MATLAB 2010b

% load swiss roll dataset with 3000 points
N = 3000;
t = (3*pi*(rand(N,1).^0.65)+pi/2);
height = 100*rand(N,1);
Y = [t.*cos(t) height  t.*sin(t)];
figure(1); scatter3(Y(:,1),Y(:,2),Y(:,3),20,t(:,1)); daspect([1 1 1]);

% compute the positive Gaussian affinities with k=20 nearest neighbors and
% sigma=5
Wp = gaussaff(Y,{'K',20},5);

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
% no more than 30 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3.
opts.maxit = 100; opts.runtime = 30; opts.tol = 1e-3;
opts.X0 = 1e-5*randn(N,d);

% Dimensionality reduction. Choose between EE, SNE or t-SNE.
% l = 1e5; [X,E,A,T] = ee(Wp,Wn,d,l,opts);
% [X,E,A,T] = ssne(Wp,d,[],opts);
[X,E,A,T] = tsne(Wp,d,[],opts);

% plot results
figure(4); scatter(X(:,1),X(:,2),24,t);
box on; daspect([1 1 1]);
title(['iter #' num2str(length(E)), ', e=' num2str(E(end)),...
       ', t=' num2str(T(end))]);

