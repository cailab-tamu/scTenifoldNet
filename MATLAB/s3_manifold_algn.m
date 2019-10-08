% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\Alignment');
% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\support');
% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\Alignment\dtw');
% addpath('C:\Users\jcai\Documents\GitHub\PCrTdMa\MATLAB\thirdparty\ManifoldWarping\Alignment\ctw');

%[aln0, aln1] = manifold_nonlinear(abs(A0),abs(A1),eye(size(A0,1)), 0.9, 10);
%return;

%%
dim=3;
mu=0.9;

W1=A0;
W2=A1;

W1=W1-diag(diag(W1));
W2=W2-diag(diag(W2));
W1=W1./max(abs(W1(:)));
W2=W2./max(abs(W2(:)));
% W1=0.5*(W1+W1');
% W2=0.5*(W2+W2');
% W1=1+W1;
% W2=1+W2;
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
vecs=vecs./vecnorm(vecs);
startx=find(vals>1e-8, 1);   % filter out eigenvalues that are ~= 0

epsilon = 1e-8;
P1 = size(W1,1);
P2 = size(W2,1);


assert(dim <= size(vecs,2)-startx+1, 'not enough eigenvectors to provide full mapping');

% Compute mappings
aln0 = vecs(1:P1,startx:dim+startx-1);
aln1 = vecs(P1+1:P1+P2,startx:dim+startx-1);
