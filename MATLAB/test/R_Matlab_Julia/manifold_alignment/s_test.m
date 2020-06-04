load X.txt
load Y.txt
[d,aln0,aln1]=i_manaln(X,Y,"R");
%[d,aln0,aln1]=i_manaln(X,Y,"Julia");


function [d,aln0,aln1]=i_manaln(A0,A1,type)
dim=30;
mu=0.9;

W1=A0+1;
W2=A1+1;

W12=eye(size(W1,2),size(W2,2));
mu = mu*(sum(W1(:))+sum(W2(:))/(2*sum(W12(:))));
W = [W1 mu*W12; mu*W12' W2];
% [~,L] = sbe_laplacian_matrix(W); to keep code consistency with python
% version, we use raw L below, instead of normalized L
% [L] = sbe_laplacian_matrix(W);
D=sum(abs(W));    %
L=diag(D)-W;

switch type
    case "R"
        [V,D] = eigs(L,dim*2,'smallestreal'); 
        d=diag(D);     % 8 min
    case "Julia"
        [V,d] = eig(L,'vector');
end            % 1.8 hrs (same Julia)  (R version 35 min)
    [d,ind] = sort(d);
    V=V(:,ind);
    V=V(:,d>=1e-8);
    V=V(:,1:dim);

    % startx=find(d>=1e-8, 1);   % filter out eigenvalues that are ~= 0
    ng=size(W1,1);
    %p2=size(W2,1);
    %assert(dim <= size(V,2)-startx+1, 'not enough eigenvectors to provide full mapping');
    %i=startx:startx+dim-1;
    %aln0 = V(1:p1,i);
    %aln1 = V(p1+1:p1+p2,i);
    aln0=V(1:ng,:);
    aln1=V(ng+1:end,:);
    d=vecnorm(aln0-aln1,2,2);
end

%[~,idx]=sort(,'descend');
% ind(1:5);
% [~,idx]=maxk(vecnorm(aln0-aln1,2,2),15);
%{
drdist=d.^2;
FC=drdist./mean(drdist);
pValues=chi2cdf(FC,1,'upper');
pAdjusted = mafdr(pValues,'BHFDR',true);
T=table([1:length(d)]',d,drdist,FC,pValues,pAdjusted);
T=sortrows(T,'d','descend');
%}
