%GO_MODULATION_BY_HOST_OF_VIRAL_PROCESS
%> A process in which a host organism modulates the frequency, rate or extent of any of a process being mediated by a virus with which it is infected. [GOC:jl]

targetg={'ANXA2','ANXA2P2','APCS','APOE','CAV2','CCL8','CCNK','CFL1','EIF2AK4','FBXL2','FMR1','IFI27','LTF','MIR155','MIR221','MIR222','NUCKS1','PC','PPIB','PRKN','PTX3','SMC3','STOM','TBC1D20','VAPA','YTHDC2','ZC3H12A','ZNF502'};
[y,clusterIdx]=ismember(targetg,genes);
clusterIdx=clusterIdx(y);

tic
m = GO.Terms(end).id;           % gets the last term id
geneschipcount = zeros(m,1);    % a vector of GO term counts for the entire chip.
genesclustercount = zeros(m,1); % a vector of GO term counts for interesting genes.

for i = 1:numel(genes)
    if isKey(SGDmap,genes{i})
        % goid = getrelatives(GO,SGDmap(genes{i}));
        % update vector counts
        geneschipcount(goid{i}) = geneschipcount(goid{i}) + 1;
        if (any(i == clusterIdx))
           genesclustercount(goid{i}) = genesclustercount(goid{i}) + 1;
        end
    end
end
if numel(unique(geneschipcount))==2, error('something wrong'); end
toc

pvalues = hygepdf(genesclustercount,max(geneschipcount),...
                  max(genesclustercount),geneschipcount);
[~,idx] = sort(pvalues);

% create a report
report = sprintf('GO Term      p-val  counts  definition\n');
for i = 1:10
    term = idx(i);
    if numel(GO(term).Terms)>0
    report = sprintf('%s%s\t%-1.4f\t%-d / %-d\t%s...\n', report, ...
                    char(num2goid(term)), pvalues(term),...
                    genesclustercount(term),geneschipcount(term),...
                    GO(term).Term.definition(2:min(end,60)));
    end
end
disp(report);



%{

topItem = idx(1);
GO(topItem).terms % the most significant gene
topItemAncestors = getancestors(GO,topItem)

secondItem = idx(2);
GO(secondItem).terms % the second most significant gene
secondItemAncestors = getancestors(GO,secondItem)

%}

%{
subGO = GO(getancestors(GO,idx(1:10)));
[cm,acc,rels] = getmatrix(subGO);
BG = biograph(cm,get(subGO.Terms,'name'));

for i = 1:numel(acc)
    pval = pvalues(acc(i));
    color = [(1-pval).^(10),pval.^(1/10),0.3];
    BG.Nodes(i).Color = color;
    BG.Nodes(i).Label = num2str(acc(i)); % add info to datatips
end

view(BG);
%}