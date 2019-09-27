function [T, W, B, Sintra, Sinter] = valid_sumsqures(data,labels,k)
%   data:    a matrix with each column representing a variable.
%   labels: a vector indicating class labels
%   W         within-group sum of squares and cross-products
%   B           between-group sum of squares and cross-products
%   T           total sum of squares and cross-products
% Sintra & Sinter: centroid diameter & linkage distance

if (size(labels, 1) == 1)
    labels = labels'; 
end

[ncase,m] = size(data);

% computing the Total sum of squares matrix
Dm = mean(data);
Dm = data - Dm(ones(ncase,1),:); 
T = Dm'*Dm;

% computing within sum of squares matrix
W = zeros(size(T));
Dm = zeros(k,m);
Sintra = zeros(1,k);
for i = 1:k
   if k > 1
      Cindex = find(labels == i);
   else
      Cindex = 1:ncase;
   end
   nk = length(Cindex);
   if nk > 1
      dataC = data(Cindex,:);
      m = mean(dataC);
      Dm(i,:) = m;
      dataC = dataC - repmat(m,nk,1);  %m(ones(nk,1),:)
      W = W + dataC'*dataC;
      dataC = sum(dataC.^2,2);
      Sintra(i) = mean(sqrt(dataC));  % distances to cluster center
   end
end

B = T - W;

% distances between cluster centers
Sinter = zeros(k,k);
if k > 1
for i = 1:k
  for j = i+1:k
     m = abs(Dm(i,:) - Dm(j,:));
     Sinter(i,j) = sqrt(sum(m.^2));
     Sinter(j,i) = Sinter(i,j);
  end
end
end