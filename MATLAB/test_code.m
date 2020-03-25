NCELLS=2000;
NGENES=400;

X0=nbinrnd(20,0.98,NGENES,NCELLS);
X0=X0(:,sum(X0)>120);
X1=X0;

X1(10,:)=X1(50,:);
X1(2,:)=X1(11,:);

genelist=strings(NGENES,1);
for k=1:NGENES, genelist(k)=sprintf("g%d",k); end

X0=sc_norm(X0,"type","libsize");
X1=sc_norm(X1,"type","libsize");

X0=zscore(X0')';
X1=zscore(X1')';

s1_network_constr;
s2_tensor_decomp;
s3_manifold_algn;

%%
d=vecnorm(aln0-aln1,2,2);
ds=d.^2;
FC=ds./mean(ds);
pValues=chi2cdf(FC,1,'upper');
pAdjusted = mafdr(pValues,'BHFDR',true);
T=table(genelist,FC,pValues,pAdjusted);

figure; 
pd=makedist('Gamma','a',0.5,'b',2);
qqplot(FC,pd);
[~,i]=sort(FC);
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1,genelist(i)};
