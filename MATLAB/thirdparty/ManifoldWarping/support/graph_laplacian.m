function L = graph_laplacian(W)
% L = graph_laplacian(W)
% W is an adjacency matrix

% L = diag(sum(W)) - W;
N = size(W,1);
D = sqrt(sum(W))';
D(D==0) = 1;
D = spdiags(1./D,0,N,N);
L = speye(N) - D*W*D;
end
