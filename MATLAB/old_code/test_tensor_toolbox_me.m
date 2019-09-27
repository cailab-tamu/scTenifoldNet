X=zeros(10,10,2,3);
for i=1:2
    for j=1:3
        X(:,:,i,j)=cov(rand(100,10));
    end
end
addpath('tensor_toolbox\');
T=tensor(X);
[M1,U1]=cp_als(T,3);
M2 = cp_als(T,3,'maxiters',100,'init',U1,'printitn',10);
fM=full(M);
