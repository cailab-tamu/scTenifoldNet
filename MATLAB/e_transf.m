function A=e_transf(A,q)
% A - adjacency matrix
if nargin<2, q=0.95; end
dim=size(A);
if numel(dim)==2
    a=max(abs(A(:)));
    if a>0
        A=A./a;
        A=A.*(abs(A)>quantile(abs(nonzeros(A)),q));
    end
elseif numel(dim)==3
    for k=1:dim(3)
        A(:,:,k)=e_transf(A(:,:,k),q);
    end
end
end