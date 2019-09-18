load pipeline_testdata.mat

addpath('ManifoldWarping\')
addpath('ManifoldWarping\Alignment\')
addpath('ManifoldWarping\support\')
addpath('ManifoldWarping\Alignment\dtw\');
addpath('ManifoldWarping\Alignment\ctw\');

X=[X0 X1];
[X,genelist]=sc_rmmtgenes(X,genelist);
[X,genelist]=sc_selectg(X,genelist,4,5);

X0=X(:,1:171);
X1=X(:,172:end);

%%

ncom=3;
tic
[A0]=sc_pcnetpar(X0,ncom);
toc
cutoff=quantile(abs(A0(:)),0.9);
A0(abs(A0)<cutoff)=0;
A0=sparse(A0);
%%
tic
[A1]=sc_pcnetpar(X1,ncom);
toc
cutoff=quantile(abs(A1(:)),0.9);
A1(abs(A1)<cutoff)=0;
A1=sparse(A1);

save pipeline_testdata.mat -append genelist2
%%
load pipeline_testdata.mat genelist2 A0 A1
clearvars -except A0 A1 genelist2
K = 8; % nearest neighbors for knn
D = 10;  % dimension of manifold

% data1 3x600
data0=A0(:,1:1500);
data1=A1(:,1:1500);

genelist=genelist2(1:1500);

[Plin,aln0,aln1,V0,V1] = manifold_warping(data0,data1,'linear',D,K,10);
