NCELLS=2000;
NGENES=400;

X0=nbinrnd(20,0.98,NGENES,NCELLS);
X0=X0(:,sum(X0)>120);
X1=X0;

X1(10,:)=X1(50,:);
X1(2,:)=X1(11,:);

%genelist=strings(NGENES,1);
%for k=1:NGENES, genelist(k)=sprintf("g%d",k); end
%sprintfc('%d',A); strsplit(num2str(A))
%
genelist=string(compose('g%d',1:NGENES)');

X0=sc_norm(X0,"type","libsize");
X1=sc_norm(X1,"type","libsize");

addpath('thirdparty\tensor_toolbox-v3.1\');
[XM0]=i_nc(X0);  % s0_network_constr;
[XM1]=i_nc(X1);  % s1_network_constr;
[A0,A1]=i_td(XM0,XM1);    % s2_tensor_decomp;

%A0=0.5*(A0+A0');
%A1=0.5*(A1+A1');

[aln0,aln1]=i_ma(A0,A1);  % s3_manifold_algn;
figure;
T=i_dr(aln0,aln1,genelist,true);   % diff regulatory gene detection

return;
%%
figure;
T2=sctenifoldnet_p(X0,X1,genelist,true);

