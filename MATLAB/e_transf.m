function A=e_transf(A)
% A - adjacency matrix
    A=A./max(abs(A(:)));
    A=A.*(abs(A)>quantile(abs(A(:)),0.95));
end