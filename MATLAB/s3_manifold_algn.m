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
W12=eye(size(W1,2),size(W2,2));

mu = mu*(sum(W1(:))+sum(W2(:))/(2*sum(W12(:))));
W = [W1 mu*W12; mu*W12' W2];

[~,L] = sbe_laplacian_matrix(W);
%%
% [vecs, vals] = eigs(L,min(dim*2,size(L,1)),'SM');
[vecs, vals] = eigs(L,2*dim,'sm');

vecs2=vecs./vecnorm(vecs);

[vals, idx] = sort(diag(vals));
vecs = vecs(:,idx);
for i=1:size(vecs,2)
    vecs(:,i) = vecs(:,i)/norm(vecs(:,i));
end


assert(isequal(vecs2,vecs))


%%
epsilon = 1e-8;
P1 = size(W1,1);
P2 = size(W2,1);

%% filter out eigenvalues that are ~= 0
for i=1:size(vals)
    if vals(i)>epsilon
        break;
    end
end
start = i;

start2=find(vals>1e-8, 1 );

assert(isequal(start, start2))


%% Compute mappings
assert(dim <= size(vecs,2)-start+1, 'not enough eigenvectors to provide full mapping');

aln0 = vecs(1:P1,start:dim+start-1);
aln1 = vecs(P1+1:P1+P2,start:dim+start-1);

% clearvars -except aln0 aln1 genelist A0 A1
