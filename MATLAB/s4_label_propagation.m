% library(apcluster)
% lk <- linKernel(MA, normalize = TRUE)
% outCluster <- apcluster(s = lk[1:100,1:100])
% save(outCluster, file = 'outCluster.RData')
% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\adapt_apV3');
% rmpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\AffinityPropagation');

[~,i]=sort(vecnorm(aln0-aln1,2,2),'descend');
gx=upper(unique(genelist(i),'stable'));

[Tpbp,Tnbp]=run_fgsea(gx,[],'bp');
[Tpmf,Tnmf]=run_fgsea(gx,[],'mf');



%%
aln=[aln0; aln1];

% s=full(compute_kernel(aln));
s=-(pdist2(aln,aln).^2);

%%
methodid=1;
switch methodid
case 1
    init_pref=median(s(s>-realmax));          % Set preference to median similarity
case 2
    N = size(aln,1);
    init_pref = min(s(s~=0));    % set the preference value as a common value of median of s
    init_pref = init_pref * log(N);
end

%%
[idx,netsim,dpsim,expref]=apcluster(s,init_pref);
C=grp2idx(idx);
% C=kmedoids(aln,700);

sortrows([unique(C) grpstats(C,C,@numel)],2)
% [length(C) sum(ans(:,2)>0) sum(ans(:,2)>1)]
sprintf('#genes=%d #modules=%d #modules (g>=2)=%d',...
    [length(C) sum(ans(:,2)>0) sum(ans(:,2)>1)]);
% clearvars -except C genelist A0 A1