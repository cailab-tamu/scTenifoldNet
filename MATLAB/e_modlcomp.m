function [G0,G1]=e_modlcomp(T,k)
if nargin<2, k=30; end

if istable(T)
    tgenes=T.genelist(T.pAdjusted<0.05);
else
    tgenes=T;
end

load(ls('A0_*.mat'),'A0','genelist');
[y]=ismember(tgenes,genelist);
assert(all(y))
[g0]=pkg.expandSeedGenesNet(A0,genelist,tgenes,k);
genelist0=genelist;

load(ls('A1_*.mat'),'A1','genelist');
[y]=ismember(tgenes,genelist);
assert(all(y))
[g1]=pkg.expandSeedGenesNet(A1,genelist,tgenes,30);
genelist1=genelist;
assert(isequal(genelist0,genelist1))

g=unique([g0;g1],'stable');
[y,i]=ismember(g,genelist);
assert(all(y))

a0=A0(i,i);
a1=A1(i,i);

G0=digraph(a0,g,'omitselfloops');
G1=digraph(a1,g,'omitselfloops');
gui.i_doublegraphs(G0,G1);
