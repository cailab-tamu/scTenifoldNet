function  [XM0,XM1]=i_nc(X0,X1,N)
if nargin<3, N=10; end

    n=size(X0,1);
    XM0=zeros(n,n,N);
    XM1=zeros(n,n,N);
    for k=1:N
        fprintf('network...%d of %d\n',k,N);

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