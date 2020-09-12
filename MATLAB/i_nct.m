function [XM]=i_nct(X,ptime,nsubsmpl,ncom,csubsmpl)
% NCT - network construction with cells subsampled using pseudotime
%
% input: X -  n (genes/features) x m (cells/samples) matrix
% input: ptime - m x 1 vector with pseudotime of cells
% output XM - k multi-layer network array (n x n x k)

if nargin<5, csubsmpl=500; end       % number of cells in subsamples                                     
if nargin<4, ncom=3; end             % number of components for PC regression
if nargin<3, nsubsmpl=10; end        % number of subsamples 


[~,t]=sort(ptime);
[n,m]=size(X);
assert(length(ptime)==m)
assert(max(t)==m)
X=X(:,t);
if m<csubsmpl*1.5, error('Too few cells.'); end
r=round(m/nsubsmpl);
winsize=max([r,csubsmpl]);
startptx=1:r:m;
while startptx(end)+winsize>m && r>1
    r=r-1;
    winsize=max([r,csubsmpl]);
    startptx=1:r:m;
    startptx=startptx(1:10);
end

    XM=zeros(n,n,nsubsmpl);
    for k=1:nsubsmpl
        fprintf('network...%d of %d\n',k,nsubsmpl);
        Xrep=X(:,startptx(k):startptx(k)+winsize);
        A=sc_pcnetpar(Xrep,ncom,true);
        XM(:,:,k)=e_transf(A);        
    end
end
