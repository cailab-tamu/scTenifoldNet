function [g1, g2] = manifold_nonlinear(W1, W2, W12, mu, dimension, epsilon)
% [g1, g2] = manifold_nonlinear(W1, W2, W12, mu, dimension, epsilon)
%
% Instance-Level Manifold Projections. Two domains.
% W1: M1*M1 matrix. weight matrix for each domain.
% W2: M2*M2 matrix
% W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
% epsilon: precision.
% mu: used to balance matching corresponding pairs and preserving manifold topology.
% dimension: dimension of input data points
% epsilon = precision. (default: 1e-8)

if nargin < 6, epsilon = 1e-8; end
P1 = size(W1,1);
P2 = size(W2,1);

%% Create weight matrix
mu = mu*(sum(sum(W1))+sum(sum(W2)))/(2*sum(sum(W12)));
W = sparse([W1 mu*W12; mu*W12' W2]);
clear W1 W2 W12;
L = graph_laplacian(W);

%% Eigen decomposition
[vecs, vals] = eigs(L,min(dimension*2,size(L,1)),'SM');
[vals, idx] = sort(diag(vals));
vecs = vecs(:,idx);
for i=1:size(vecs,2)
    vecs(:,i) = vecs(:,i)/norm(vecs(:,i));
end

%% filter out eigenvalues that are ~= 0
for i=1:size(vals);
    if vals(i)>epsilon
        break;
    end
end
start = i;

%% Compute mappings
assert(dimension <= size(vecs,2)-start+1, 'not enough eigenvectors to provide full mapping');

g1 = vecs(1:P1,       start:dimension+start-1);
g2 = vecs(P1+1:P1+P2, start:dimension+start-1);

end


