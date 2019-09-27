function [M, dmax, Dist] = simatrix_make(data,type,nrow)
%   data:    a matrix with each column representing a variable.

if type == 1
    [Dist, dmax] = similarity_euclid(data,2);   %  pdist(data,'type');
else
    Dist = 1-(1+similarity_pearson(data'))/2;
    dmax = 1;
end

nap = nrow*nrow-nrow;
M = zeros(nap,3);
j = 1;
for i=1:nrow
    for k = [1:i-1,i+1:nrow]
        M(j,1) = i;
        M(j,2) = k;
        M(j,3) = -Dist(i,k);
        j = j+1;
    end
end