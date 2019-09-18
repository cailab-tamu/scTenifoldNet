function idx = knnsearch_same(X,k)
% same as idx = knnsearch(X,[],k), but faster
aa = sum(X.*X);
n = size(aa,2);
D = repmat(aa',[1 n]) + repmat(aa,[n 1]) - 2*X'*X;
[~,idx] = sort(D);
idx = idx(2:k+1,:)';
