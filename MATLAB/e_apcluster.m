function [C,GS]=e_apcluster(aln0,genelist)

if nargin<2, genelist=[]; end

% s=full(compute_kernel(aln));
s=-(pdist2(aln0,aln0).^2);

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
a=sortrows([unique(C) grpstats(C,C,@numel)],2,'descend');
% [length(C) sum(ans(:,2)>0) sum(ans(:,2)>1)]
fprintf('#genes=%d #modules=%d #modules (g>=2)=%d',...
    [length(C) sum(a(:,2)>0) sum(a(:,2)>1)]);
GS=cell(size(a,1),1);
if ~isempty(genelist)
for k=1:size(a,1)
    GS{k}=genelist(C==a(k,1));
    % GS{k}=sprintf('%s,',genelist(C==a(k,1)));
end
end