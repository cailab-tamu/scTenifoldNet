gunzip('goa_human.gaf.gz')
T=readtable('goa_human.gaf','filetype','text');
DB_Object_Symbol=T.Var3;
GOid=T.Var5;
delete('goa_human.gaf')

tic
SGDmap = containers.Map();
for i = 1:length(DB_Object_Symbol)
    key = DB_Object_Symbol{i};
    a=str2num(GOid{i}(4:end));
    if isKey(SGDmap,key)
        SGDmap(key) = [SGDmap(key) a];
    else
        SGDmap(key) = a;
    end
end
toc

fprintf('Number of annotated genes related to molecular function is %d.\n',SGDmap.Count)
fprintf('Number of unique GO terms associated to annotated genes is %d.\n',numel(unique(GOid)))
fprintf('Number of gene-GO term associations is %d.\n',numel(DB_Object_Symbol))

GO = geneont('live',true); % this step takes a while
% get(GO)

%%
% clusterIdx=1:50;
genes=unique(DB_Object_Symbol);

tic
for i = 1:numel(genes)
    if isKey(SGDmap,genes{i})
        goid{i} = getrelatives(GO,SGDmap(genes{i}));
    end
end
toc

