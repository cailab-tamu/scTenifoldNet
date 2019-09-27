load neuron_microglia.mat Xne agegroup_ne genelist

c=grp2idx(agegroup_ne);
X0=Xne(:,c==1);
X1=Xne(:,c==2);

i0=mean(X0>0,2)>0.1;
i1=mean(X1>0,2)>0.1;
genelist0=genelist(i0);
genelist1=genelist(i1);
X0=X0(i0,:);
X1=X1(i1,:);

[genelist,i,j]=intersect(genelist0,genelist1);
X0=X0(i,:);
X1=X1(j,:);

% clearvars -except X0 X1 agegroup_ne genelist




%%
for k=1:10
    Xrep=X0(:,randperm(size(X0,2)));    
    A=sc_pcnetpar(Xrep(:,1:700),3,true);
    Ayng{k}=sparse(A.*(abs(A)>quantile(abs(A(:)),0.95)));
end
save neuron_tensor Ayng genelist

%%
for k=1:10
    k
    tic
    Xrep=X1(:,randperm(size(X1,2)));
    A=sc_pcnetpar(Xrep(:,1:700),3,true);
    Aold{k}=sparse(A.*(abs(A)>quantile(abs(A(:)),0.95)));
    toc
end

save neuron_tensor -append Aold





