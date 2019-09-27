function [clusters, newlabels, Clength] = ind2cluster(labels)

C = unique(labels);
newlabels = labels;
k = length(C);
clusters = cell(1,k);
Clength = zeros(1,k);

for i = 1:k
  ind = find(labels==C(i));
  clusters{i} = ind;
  newlabels(ind) = i;
  Clength(i) = length(ind);
end