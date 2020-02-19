% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\Alignment');
% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\support');
% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\Alignment\dtw');
% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\Alignment\ctw');
%[aln0, aln1] = manifold_nonlinear(abs(A0),abs(A1),eye(size(A0,1)), 0.9, 10);
%return;

%%
dim=30;
mu=0.9;

W1=A0+1;
W2=A1+1;

% W1=W1-diag(diag(W1));
% W2=W2-diag(diag(W2));
%W1=W1./max(abs(W1(:)));
%W2=W2./max(abs(W2(:)));

% W1(W1~=0)=W1(W1~=0)+1;
% W2(W2~=0)=W2(W2~=0)+1;
% W1(W1<0)=0;
% W2(W2<0)=0;

W12=eye(size(W1,2),size(W2,2));
mu = mu*(sum(W1(:))+sum(W2(:))/(2*sum(W12(:))));

W = [W1 mu*W12; mu*W12' W2];

% [~,L] = sbe_laplacian_matrix(W); to keep code consistency with python
% version, we use raw L below, instead of normalized L
[L] = sbe_laplacian_matrix(W);

%%
% [vecs, vals] = eigs(L,min(dim*2,size(L,1)),'SM');

%{
errortag=false;
lastwarn('');
try
    [vecs, vals] = eigs(L,2*dim,'sm');
catch    
    errortag=true;
end
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg), errortag=true; end
if errortag
    disp('Now using eig.m...');
    [vecs, vals] = eig(L);
    [vals, idx] = sort(diag(vals));
    vecs = vecs(:,idx);
end
%}

% [V,D]=eigs(L,dim*2,'SM');
[V,D] = eig(L);
[d,ind] = sort(diag(D));
V=V(:,ind);
V=V(:,d>=1e-8);
V=V(:,1:dim);

% startx=find(d>=1e-8, 1);   % filter out eigenvalues that are ~= 0
p1=size(W1,1);
% p2=size(W2,1);
%assert(dim <= size(V,2)-startx+1, 'not enough eigenvectors to provide full mapping');
%i=startx:startx+dim-1;
%aln0 = V(1:p1,i);
%aln1 = V(p1+1:p1+p2,i);
aln0=V(1:p1,:);
aln1=V(p1+1:end,:);
