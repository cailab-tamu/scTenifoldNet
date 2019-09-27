function Outs = valid_external(index1,c2)

Outs = [];
N = size(index1,2);
Q = sum(double(int32(c2))-c2) ~= 0;
Q = Q || sum(sum(double(int32(index1))-index1)) ~= 0;
if Q
   return;
end

for i = 1:N
c1 = index1(:,i);
C=Contingency(c1,c2);	  %form contingency matrix
n = length(c1);
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

% nchoosek(n,2); combntns(1:5,3);
ns = n*(n-1)/2;	                    %total number of pairs        
sumC=sum(sum(C.^2));	 %sum over rows & columnns of nij^2
sumij = nis+njs;

R = ns+sumC-sumij*0.5;    %no. agreements; disagreements -t2+sumij*0.5;

nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));
if ns==nc
   AR=0;			                   % if k=1, define Rand = 0
else
   AR=(R-nc)/(ns-nc);	       % adjusted Rand - Hubert & Arabie 1985
end

Rand = R/ns;                       % Rand 1971
Jac = (sumC-n)/(sumij-sumC-n);	 % Jaccard - Jain and Dubes 1988

ni = sum(C,2);
ni = ni.*(ni-1)/2;
nis = sum(ni);
nj = sum(C,1);
nj = nj.*(nj-1)/2;
njs = sum(nj);
FM = 0.5*(sumC-n)/sqrt(nis*njs);    % FM - Fowlkes and Mallows 1983

Outs = [Outs [Rand; AR; Jac; FM]];
end

function Cont = Contingency(Mem1,Mem2)
Cont = zeros(max(Mem1),max(Mem2));
for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i)) = Cont(Mem1(i),Mem2(i))+1;
end
