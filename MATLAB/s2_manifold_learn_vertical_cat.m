addpath('ManifoldWarping\')
addpath('ManifoldWarping\Alignment\')
addpath('ManifoldWarping\support\')
addpath('ManifoldWarping\Alignment\dtw\');
addpath('ManifoldWarping\Alignment\ctw\');

load pipeline_testdata.mat genelist2 A0 A1 X0 X1
% clearvars -except A0 A1 genelist2
K = 8; % nearest neighbors for knn
D = 15;  % dimension of manifold

n=11697;

% data1 3x600
W0=A0(1:n,1:n);
W1=A1(1:n,1:n);
% X0=double(X0);
% X1=double(X1);
% X0=run_magic(X0);
% X1=run_magic(X1);

X0=double(X0(1:n,:))';
X1=double(X1(1:n,:))';

genelist=genelist2(1:n);

% [~,aln0,aln1] = manifold_warping_given_W(X0,X1,W0,W1,'linear',D,K,10);
% [~,aln0,aln1] = manifold_warping(W0,W1,'linear',D,K,10);

[~,aln0,aln1] = manifold_warping(X0,X1,'linear',D,K,10);


%%
aln=[aln0'; aln1'];

C=kmedoids(aln,700);
clc
gg=[genelist;genelist];
for k=1:max(C)
    i=C==k;
    fprintf('%d %d ',k,sum(i));
    fprintf('%s ',gg(i))
    fprintf('\n');
end
%%
P_v=[];
G0_v=[]; G1_v=[];
J_v=[];

for k=1:max(C)
    
idx=find(C==k);
if length(idx)>20
    gid=idx;
    g0id=gid(gid<=n);
    g1id=gid(gid>n)-n;
    jacx=length(intersect(g0id,g1id))/length(union(g0id,g1id));
    if length(g0id)>length(g1id)
        J_v=[J_v; jacx];
    else
        J_v=[J_v; -jacx];
    end
    if ~isempty(g0id) && ~isempty(g1id)
        B0 = reshape(A0(g0id,g0id),[],1);
        B1 = reshape(A1(g1id,g1id),[],1);    
        [~,p]=kstest2(full(B0),full(B1));        
        P_v=[P_v;p];
    else
        P_v=[P_v;nan];
    end
    G0_v=[G0_v;string(sprintf('%s ',genelist(g0id)))];
    G1_v=[G1_v;string(sprintf('%s ',genelist(g1id)))];
end
end

absJ_v=abs(J_v);
Tres=table(absJ_v,J_v,P_v,G0_v,G1_v);
writetable(Tres,'res_v','filetype','spreadsheet');



