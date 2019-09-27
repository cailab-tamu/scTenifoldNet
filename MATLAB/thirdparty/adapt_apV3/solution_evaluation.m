function [NC,Sil,Silmin,NCfix] = solution_evaluation(data,M,labels,NC,NCfix,simatrix,nrow,type,cut)

if type == 1
    type = 'euclidean';
else
    type = 'correlation';
end

dim = 0;
if simatrix
%   M = simatrix_make(data,type,nrow);
    Ms = NaN*ones(nrow,nrow);
    for i=1:nrow
        Ms(i,i) = 0;
    end
    for i=1:size(M,1)
        ni=M(i,1);
        nj=M(i,2);
        Y = -M(i,3);
        Ms(ni,nj)= Y;
        if Y < dim
            dim = Y;
        end
    end
end

if dim < 0
    Ms = Ms - dim + 1;
    for i=1:nrow
        Ms(i,i) = 0;
    end
end

dim = length(NC);
Sil =[];
Sildelete = [];
for i = 1:dim
    Y = labels(:,i);
    if simatrix
        Smax = silhouette2(data, Y, Ms);
    else
        Smax = silhouette(data, Y, type);
    end
    dn = isfinite(Smax);
    Sil(i) = mean(Smax(dn));
    [C, Y, dmax] = ind2cluster(Y);
    dmax = min(dmax);
    Sildelete(i) = dmax < cut;
    Q =[];
    for j = 1:length(C)
        R = C{j};
        R = Smax(R);
        dn = isfinite(R);
        Q(j) = mean(R(dn));
    end
    Silmin(i) = min(Q);
end

Sildelete = find(Sildelete);
if length(Sildelete) < dim
    Sil(Sildelete) = [];
    Silmin(Sildelete) = [];
    NC(Sildelete) = [];
    NCfix(Sildelete) = [];
end
