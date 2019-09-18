function A=sbe_affinity_matrix(X,methodid)
if nargin<2
    methodid=1;
end
switch methodid
    case 1
        A=full(compute_alpha_kernel_sparse(X));
    case 2
        A=full(compute_kernel_sparse(X));
    case 3
        D=pdist2(X,X);
        A=scale_dist3(D,7);
    case 4
        D=pdist2(X,X);
        kn=floor(log2(size(X,1)))+1;
        A=scale_dist3_knn(D,7,kn,false);
end

