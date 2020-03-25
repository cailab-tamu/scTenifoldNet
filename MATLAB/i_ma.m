function [aln0,aln1]=i_ma(A0,A1,dim)

if nargin<3, dim=30; end
mu=0.9;

W1=A0+1;
W2=A1+1;

W12=eye(size(W1,2),size(W2,2));
mu = mu*(sum(W1(:))+sum(W2(:))/(2*sum(W12(:))));
W = [W1 mu*W12; mu*W12' W2];

D=sum(abs(W));
L=diag(D)-W;

[V,D] = eigs(L,dim*2,'smallestreal'); d=diag(D);
% [V,d] = eig(L,'vector');
[d,ind] = sort(d);
V=V(:,ind);
V=V(:,d>=1e-8);
V=V(:,1:dim);

p1=size(W1,1);
aln0=V(1:p1,:);
aln1=V(p1+1:end,:);
