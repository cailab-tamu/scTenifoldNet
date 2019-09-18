X=X';
kneighbors=5;
[dim, N] = size(X);
[kNNgraphlength, Graph] = kNNgraphmex(X(:), N, dim, kneighbors, 1);
C = reshape(Graph, kneighbors+1, N );
