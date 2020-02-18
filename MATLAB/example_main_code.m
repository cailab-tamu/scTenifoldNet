addpath('../../scGEAToolbox/');
load example10xdata.mat
[X,genelist]=sc_selectg(X,genelist,10,7);
[X]=sc_selectc(X,2000);

% rng default
c=rand(size(X,2),1)>0.5;
X0=X(:,~c);
X1=X(:,c);
X0=X1;
i=find(genelist=="SWAP70");    % 1640
j=find(genelist=="CISD1");     % 301
s=X1(i,:); t=X1(j,:);
X1(i,:)=t; X1(j,:)=s;
%%
s1_network_constr;
s2_tensor_decomp;
s3_manifold_algn;

[~,i]=sort(vecnorm(aln0-aln1,2,2),'descend');
gx=upper(unique(genelist(i),'stable'));

find(gx=="SWAP70")
find(gx=="CISD1")

% s4_label_propagation;
% s5_module_gene_anno;
