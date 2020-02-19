%addpath('../../scGEAToolbox/example_data/example10xdata.mat');
%[X,genelist]=sc_selectg(X,genelist,10,7);
%[X]=sc_selectc(X,2000);
load example10xdata_1877g X genelist
%%
rng default
c=rand(size(X,2),1)>0.5;
X0=X(:,~c);
X1=X(:,c);
X0=X1;
i=genelist=="SWAP70";    % 8
j=genelist=="CISD1";     % 80
% s=X1(i,:); t=X1(j,:);
% X1(i,:)=t; X1(j,:)=s;
X1(i,:)=X1(j,:);
%%
s1_network_constr;
s2_tensor_decomp;
s3_manifold_algn;
[~,i]=sort(vecnorm(aln0-aln1,2,2),'descend');
gx=genelist(i);
find(gx=="SWAP70")
find(gx=="CISD1")
