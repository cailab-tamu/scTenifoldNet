function [map1, map2] = manifold_linear(X1, X2, W1, W2, W12, mu, max_dim, epsilon)
% [map1, map2] =  manifold_linear(X1, X2, W1, W2, W12, [mu, max_dim, epsilon])
%
% Feature-level Manifold Projections. Two domains.
% X1: P1*M1 matrix, M1 examples in a P1 dimensional space.
% X2: P2*M2 matrix
% W1: M1*M1 matrix. weight matrix for each domain.
% W2: M2*M2 matrix
% W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
% mu: used to balance matching corresponding pairs and preserving manifold topology.
% max_dim: max dimensionality of the new space. (default: 200)
% epsilon: precision. (default: 1e-8)

if nargin < 8, epsilon = 1e-8; end
if nargin < 7, max_dim = 200; end
if nargin < 6, mu = 1; end

%% get sizes for convenience later
[P1,M1] = size(X1);
[P2,M2] = size(X2);
%% Create weight matrix
mu = mu*(sum(sum(W1))+sum(sum(W2)))/(2*sum(sum(W12)));
W = sparse([W1 mu*W12; mu*W12' W2]);
clear W1 W2 W12;
L = graph_laplacian(W);

%% prepare for decomposition
Z = sparse([X1 zeros(P1,M2); zeros(P2, M1) X2]);
[u, s, ~] = svds(Z*Z',max_dim+1);
Fplus = pinv(u*sqrt(s));
T = Fplus*Z*L*Z'*Fplus';

clear u s Z W L;

%% Eigen decomposition
[vecs, vals] = eigs((T+T')/2,min(max_dim,size(T,1)-1),'SM');
[vals, idx] = sort(diag(vals));
vecs = Fplus'*vecs(:,idx);
clear T Fplus;
for i=1:size(vecs,2)
    vecs(:,i) = vecs(:,i)/norm(vecs(:,i));
end

%% filter out eigenvalues that are ~= 0
for i=1:size(vals)
    if vals(i)>epsilon
        break;
    end
end
start = i;

%% Compute mappings
m = min(max_dim,size(vecs,2)-start+1);
map1 = vecs(1:P1,       start:m+start-1);
map2 = vecs(P1+1:P1+P2, start:m+start-1);

end


