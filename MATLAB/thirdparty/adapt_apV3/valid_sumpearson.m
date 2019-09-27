function [T, W, B, Sintra, Sinter] = valid_sumpearson(data,labels,k)
% within-, between-cluster and total sum of squares
% Sintra/Sinter: centroid diameter/linkage distance based on Pearson correlation

C = mean(data);
R = similarity_pearsonC(data', C');
T = R*R';

W = 0;
Sintra = zeros(1,k);
Sinter = zeros(k,k);
for i = 1:k
   Ui = find(labels == i);
   ni = length(Ui);
   if ni > 1
      datai = data(Ui,:);
      C = mean(datai);
      R = similarity_pearsonC(datai', C'); % distances to cluster center
      Sintra(i) = mean(R); 
      W = W + R*R';
   end
   % distances between cluster centers
  for j = i+1:k
     Ui = find(labels == j);
     ni = length(Ui);
     if ni > 0
        datai = data(Ui,:);
        if ni == 1
          Cj = datai;
        else
          Cj = mean(datai);
        end
        Sinter(i,j) = similarity_pearsonC(Cj', C');
     end
  end
end

B = T - W;