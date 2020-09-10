function [A]=sc_pcnetpar(X,ncom,fastersvd)
% [A]=sc_pcnetpar(X,ncom)
% ncom - number of components used (default=3)
% ref: https://rdrr.io/cran/dna/man/PCnet.html
% https://github.com/cran/dna/blob/master/src/rpcnet.c
% https://rdrr.io/cran/dna/f/inst/doc/Introduction.pdf

if nargin<2, ncom=3; end
if nargin<3, fastersvd=false; end
opts.maxit=150;

if fastersvd
    opts.maxit=150;
    pw1=fileparts(which(mfilename));
    pth=fullfile(pw1,'thirdparty/faster_svd/lmsvd');
    addpath(pth);
end

% [X]=sc_norm(X);
X=X';
X=zscore(X);

n=size(X,2);
A=1-eye(n);
B=A(:,1:end-1);

parfor k=1:n
    y=X(:,k);
    Xi=X;
    Xi(:,k)=[];    
    if fastersvd
        warning off
        [~,~,coeff]=lmsvd(Xi,ncom,opts);
        warning on
    else
        [~,~,coeff]=svds(Xi,ncom);
    end    
    score=Xi*coeff;
    score=(score./(vecnorm(score).^2));
    Beta=sum(y.*score);
    B(k,:)=coeff*Beta';
end
for k=1:n
   A(k,A(k,:)==1)=B(k,:); 
end
