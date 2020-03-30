function T=sctenifoldnet_m(X0,X1,genelist,doplot)
    if nargin<3
        error('USAGE: T=sctenifoldnet_m(X0,X1,genelist)');
    end
    if nargin<4, doplot=false; end
    if size(X0,1)~=size(X1,1)
        error('X0 and X1 need the same number of rows.');
    end
    if size(X0,1)~=length(genelist)
        error('Length of genelist should be the same as the number of rows of X0 or X1.');
    end
    pw0=pwd;
    pw1=fileparts(which(mfilename));
    cd(pw1);
    addpath('thirdparty\tensor_toolbox-v3.1\');
    
    n=size(X0,1);
    XM0=zeros(n,n,10);
    XM1=zeros(n,n,10);
    for k=1:10
        fprintf('network...%d of 10\n',k);

        Xrep=X0(:,randperm(size(X0,2)));
        A=sc_pcnetpar(Xrep(:,1:500),3,false);
        A=A./max(abs(A(:)));        
        XM0(:,:,k)=A.*(abs(A)>quantile(abs(A(:)),0.95));

        Xrep=X1(:,randperm(size(X1,2)));
        A=sc_pcnetpar(Xrep(:,1:500),3,false);    
        A=A./max(abs(A(:)));
        XM1(:,:,k)=A.*(abs(A)>quantile(abs(A(:)),0.95));
    end
        
    A0=i_td(XM0,2);
    A1=i_td(XM1,2);
    
    [aln0,aln1]=i_ma(A0,A1);
    drdist=vecnorm(aln0-aln1,2,2).^2;
    FC=drdist./mean(drdist);
    pValues=chi2cdf(FC,1,'upper');
    pAdjusted = mafdr(pValues,'BHFDR',true);
    if size(genelist,1)==1, genelist=genelist'; end
    T=table(genelist,drdist,FC,pValues,pAdjusted);    
    if doplot
        pd = makedist('Gamma','a',0.5,'b',2);
        qqplot(FC,pd);
        [~,i]=sort(FC);
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist(i)};
    end
    cd(pw0);
end

