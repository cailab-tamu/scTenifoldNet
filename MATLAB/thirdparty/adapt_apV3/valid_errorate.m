function Rerror = valid_errorate(labels, truelabels)
% computing error rates for every clusters if true labels are given

nrow = length(truelabels);
[R, truelabels] = ind2cluster(truelabels);
[truelabels, S] = sort(truelabels);
labels = labels(S);
R = ind2cluster(labels);
k = length(R);

S = ones(nrow,1);
Sm=[];
for i = 1:k
   Sk = sort(R{i});
   m = length(Sk);
   n = max([round(sqrt(sqrt(m))+0.1) 1]);
   Sm(i) = mean(Sk(n:m-n+1));
end
[Sk, Sm]=sort(Sm);

for i = 1:k
   Sk = Sm(i);
   Sk = R{Sk};
   S(Sk)=i;
   Q = S(Sk)-truelabels(Sk);
   Q = nonzeros(Q);
   Q = length(Q);
   Rerror = 100*Q/length(Sk);
   fprintf('\n Error rate of cluster %d : %4.2f %%',i, Rerror);
end

S=S-truelabels;
S=nonzeros(S);
S=length(S);
Rerror=100*S/nrow;
fprintf('\n Error rate for all the data: %4.2f %% \n',Rerror);
