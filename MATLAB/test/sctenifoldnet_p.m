function T=sctenifoldnet_p(X0,X1,genelist,doplot)
    if nargin<3
        error('USAGE: T=sctenifoldnet_p(X0,X1,genelist)');
    end
    if nargin<4, doplot=false; end
    if size(X0,1)~=size(X1,1)
        error('X0 and X1 need the same number of rows.');
    end
    if size(X0,1)~=length(genelist)
        error('Length of genelist should be the same as the number of rows of X0 or X1.');
    end    
    a0=sc_pcnetpar(X0,3,false);
    a1=sc_pcnetpar(X1,3,false);
    if exist('svdsketch.m','file')
        [U,S,V]=svdsketch(a0);
    else
        [U,S,V]=svds(a0,50);
    end
    A0=U*S*V';
    
    if exist('svdsketch.m','file')
        [U,S,V]=svdsketch(a1);
    else
        [U,S,V]=svds(a1,50);
    end    
    A1=U*S*V';
    [aln0,aln1]=i_ma(A0,A1);
    drdist=vecnorm(aln0-aln1,2,2);
    %drdist=drdist.^2;
    FC=drdist./mean(drdist);
    pValues=chi2cdf(FC,1,'upper');
    pAdjusted = mafdr(pValues,'BHFDR',true);
    if size(genelist,1)==1, genelist=genelist'; end
    
   sortid=(1:length(genelist))';
   if size(genelist,2)>1, genelist=genelist'; end
   T=table(sortid,genelist,drdist,FC,pValues,pAdjusted);
   T = sortrows(T,'drdist','descend');

   if doplot
        pd = makedist('Gamma','a',0.5,'b',2);
        qqplot(FC,pd);
        [~,i]=sort(FC);
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist(i)};
    end
end

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
end

function txt = i_myupdatefcn1(~,event_obj,g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
end