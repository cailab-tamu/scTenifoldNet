function  [XM0,XM1]=i_nc(X0,X1)
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
end