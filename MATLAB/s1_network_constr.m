% load neuron_microglia.mat Xne agegroup_ne genelist
% Xinput=Xne;
% agegroup_input=agegroup_ne;
% c=grp2idx(agegroup_input);


%{
i0=mean(X0>0,2)>0.1;
i1=mean(X1>0,2)>0.1;
genelist0=genelist(i0);
genelist1=genelist(i1);
X0=X0(i0,:);
X1=X1(i1,:);

[genelist,i,j]=intersect(genelist0,genelist1);
X0=X0(i,:);
X1=X1(j,:);
%}

% clearvars -except X0 X1 agegroup_ne genelist

%%
for k=1:10
    k
    tic    
    Xrep=X0(:,randperm(size(X0,2)));    
    A=sc_pcnetpar(Xrep(:,1:500),3,true);
    A0{k}=sparse(A.*(abs(A)>quantile(abs(A(:)),0.95)));
    toc
end

%%
for k=1:10
    k
    tic
    Xrep=X1(:,randperm(size(X1,2)));
    A=sc_pcnetpar(Xrep(:,1:500),3,true);
    A1{k}=sparse(A.*(abs(A)>quantile(abs(A(:)),0.95)));
    toc
end






