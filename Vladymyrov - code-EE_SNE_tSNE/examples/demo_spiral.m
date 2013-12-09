% Example file for embedding 2D Spiral dataset into 1D

% we fix a random seed for repeatability of the results
s = RandStream('mcg16807','Seed',29); %RandStream.setGlobalStream(s);
RandStream.setDefaultStream(s);

% generate spiral dataset with 200 points
N = 200;
a = pi/2; b = 4*pi;
t = ((b-a)*(linspace(0,1,N).^0.65)+a)';
Y = [t.*sin(t) t.*cos(t)];
figure(1); plot(Y(:,1),Y(:,2),'o');
box on; daspect([1 1 1]); 
title('Spiral dataset')

% compute the positive Gaussian affinities with k=5 nearest neighbors and
% sigma=1
Wp = gaussaff(Y,{'K',5},1);

% compute the negative Gaussian affinities (for EE only) based on the
% square distance matrix
Wn = sqdist(Y);

% Make sure that both matrices are normalized, symmetric and have zeros
% on the diagonal
Wp = (Wp+Wp')/2; Wp(1:N+1:N^2) = 0; Wp = Wp/sum(Wp(:));
Wn = (Wn+Wn')/2; Wn(1:N+1:N^2) = 0; Wn = Wn/sum(Wn(:));

% we want one dimensional embedding
d = 1;

% stopping criteria: number of iterations is no more than 100, runtime is
% no more than 5 seconds, and the relative tolerance in the embedding is 
% no less than 1e-3.
opts.maxit = 100; opts.runtime = 5; opts.tol = 1e-3;

% Dimensionality reduction. Choose between EE, SNE or t-SNE.
l = 1; [X,E,A,T] = ee(Wp,Wn,d,l,opts);
% [X,E,A,T] = ssne(Wp,d,[],opts);
% [X,E,A,T] = tsne(Wp,d,[],opts);

% plot results
figure(3); plot(1:N,X,'bo');
title(['iter #' num2str(length(E)), ', e=' num2str(E(end)),...
       ', t=' num2str(T(end))]);

