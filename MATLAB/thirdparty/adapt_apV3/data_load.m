p = [];
if simatrix == 1
    Ms = load(sw);
    if id == 12 || id == 13
        M = Ms.s;       % given similarity matrix
        p = Ms.p;
        name = Ms.x;
    elseif id == 14
        M = Ms.s;
        p = [];
        refseq_exon = Ms.refseq_exon;
        refseq_intron = Ms.refseq_intron;
        [M,dim] = matrix_transform(M,nsubset,nrow);
        refseq_exon = refseq_exon(1:nsubset-1);
        refseq_intron = refseq_intron(1:nsubset-1);
        nrow = nsubset;
    else
        M = Ms;
    end
    Ms = [];
    data = [];
else
    data = load(sw);
    [nrow, dim] = size(data);
    M = [];
end

% taking true class labels from a data file
truelabels = ones(nrow,1);
if id < 11 || (id > 20 && id < 30)  % when 1st column is class labels
   truelabels = data(:,1); 
   data = data(:,2:dim);
   dim = dim-1;
end
