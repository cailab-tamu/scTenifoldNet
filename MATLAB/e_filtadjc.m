function A=e_filtadjc(A,q)
% A - adjacency matrix
if nargin<2, q=0.95; end
dim=size(A);
if numel(dim)==2
    a=mean(maxk(abs(A(:)),10));     % top 10 average
    if a>0
        A=A./a;
        A=A.*(abs(A)>quantile(abs(nonzeros(A)),q));
    end
elseif numel(dim)==3
    for k=1:dim(3)
        A(:,:,k)=e_filtadjc(A(:,:,k),q);
    end
end
end
