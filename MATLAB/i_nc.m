function  [XM0,XM1]=i_nc(X0,X1,N,ncom,usingboot)
if nargin<5, usingboot=false; end   % using m-out-of-n bootstrap (false by default)
if nargin<4, ncom=3; end    % number of components
if nargin<3, N=10; end      % number of subsamples 

    n=size(X0,1);
    XM0=zeros(n,n,N);
    XM1=zeros(n,n,N);
    for k=1:N
        fprintf('network...%d of %d\n',k,N);
        
        n0=size(X0,2);
        Xrep=X0(:,randperm(n0));
        if usingboot, i=randi(n0,1,500); Xrep=Xrep(:,i); end
        A=sc_pcnetpar(Xrep(:,1:500),ncom,false);
        A=A./max(abs(A(:)));        
        XM0(:,:,k)=A.*(abs(A)>quantile(abs(A(:)),0.95));

        n1=size(X1,2);
        Xrep=X1(:,randperm(n1));
        if usingboot, i=randi(n1,1,500); Xrep=Xrep(:,i); end
        A=sc_pcnetpar(Xrep(:,1:500),ncom,false);
        A=A./max(abs(A(:)));
        XM1(:,:,k)=A.*(abs(A)>quantile(abs(A(:)),0.95));
    end
end