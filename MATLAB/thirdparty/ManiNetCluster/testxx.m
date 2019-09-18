% T1=readtable('dayOrthoExpr.csv');
% T2=readtable('nightOrthoExpr.csv');
% X=table2array(T1(:,2:end));
% Y=table2array(T2(:,2:end));
% genelist=string(T1.Var1);
% assert(isequal(genelist,string(T2.Var1)))
% save maninetcluster_inputdata X Y genelist

% reticulate::use_python("C:\\Program Files\\Python36")
% reticulate::use_python("C:\\Users\\jcai\\AppData\\Local\\Continuum\\miniconda3")
% reticulate::py_config()
load maninetcluster_inputdata X Y genelist
% A=full(compute_alpha_kernel_sparse(X));
addpath('ManifoldWarping\matlab\Alignment\');
addpath('ManifoldWarping\matlab\Alignment\dtw\');
addpath('ManifoldWarping\matlab\Alignment\ctw\');
addpath('ManifoldWarping\matlab\support\');

%%
Wx=sbe_affinity_matrix(X);
Wy=sbe_affinity_matrix(Y);
%Dx=diag(sum(Wx));
%Dy=diag(sum(Wy));
W=eye(size(Wx,1));
Z=blkdiag(X,Y);
% D=blkdiag(Dx,Dy);
W=0.5*[Wx, W; W' Wy];
% L=D-W;
[~,L]=sbe_laplacian_matrix(W);
max_dim=200;


%  Z = bdiag(X,Y)
%   mu = (sum(Wx)+sum(Wy))/(2*sum(sum(W)));
%   
%   Wall = rbind(cbind(Wx,W*mu),cbind(t(W*mu),Wy))
%   Ddiag = rowSums(Wall)
%   D = diag(rowSums(Wall), nrow = mx+my, ncol = mx+my)
%   L = D - Wall
%   
%   A1 = as.matrix(t(Z) %*% L %*% Z)
%   A2 = t(Z) %*% D  %*% Z
%   # A2 = (A2 + t(A2))/2
%   ## here is a dirty quick solution, it need that both A1 positive def
%   MP <- solve(A1, A2)
%   MP <- eigs(MP, k = d)
%   P = Re(MP$vectors)



[u, s, ~] = svds(Z*Z',max_dim+1);
Fplus = pinv(u*sqrt(s));
T = Fplus*Z*L*Z'*Fplus';

%% Eigen decomposition
[vecs, vals] = eigs((T+T')/2,min(max_dim,size(T,1)-1),'SM');
[vals, idx] = sort(diag(vals));
vecs = Fplus'*vecs(:,idx);
% clear T Fplus;
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


%%
K = 8; % nearest neighbors for knn
D = 3;  % dimension of manifold
[Pnon,aln1,aln2] = manifold_warping(X',Y','linear',D,K,10);






