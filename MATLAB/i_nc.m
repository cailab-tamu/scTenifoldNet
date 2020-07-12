function  [XM0,XM1]=i_nc(X0,X1,N,ncom,usebootstrp)

if nargin<5, usebootstrp=false; end   % using m-out-of-n bootstrap (false by default)
                                    % using jackknife (by default)
if nargin<4, ncom=3; end    % number of components
if nargin<3, N=10; end      % number of subsamples 

    n=size(X0,1);
    XM0=zeros(n,n,N);
    XM1=zeros(n,n,N);
    for k=1:N
        fprintf('network...%d of %d\n',k,N);
        
        n0=size(X0,2);
        if usebootstrp % bootstrap 
            i=randi(n0,1,500);
            Xrep=X0(:,i);
        else     % jackknife
            Xrep=X0(:,randperm(n0));
            Xrep=Xrep(:,1:500);
        end    
        
        %if any(sum(Xrep,1)==0)||any(sum(Xrep,2)==0)
        %    Xrep=i_goodrep(X0,usingboot);
        %end
        
        A=sc_pcnetpar(Xrep,ncom,false);
        A=A./max(abs(A(:)));        
        XM0(:,:,k)=A.*(abs(A)>quantile(abs(A(:)),0.95));

        n1=size(X1,2);        
        if usebootstrp % bootstrap 
            i=randi(n1,1,500);
            Xrep=X1(:,i);
        else     % jackknife
            Xrep=X1(:,randperm(n1));
            Xrep=Xrep(:,1:500);
        end
        %if any(sum(Xrep,1)==0)||any(sum(Xrep,2)==0)
        %    Xrep=i_goodrep(X1,usingboot);
        %end
        A=sc_pcnetpar(Xrep,ncom,false);
        A=A./max(abs(A(:)));
        XM1(:,:,k)=A.*(abs(A)>quantile(abs(A(:)),0.95));
    end
end

function Xrep=i_goodrep(X,usingboot)
    disp('using i_goodrep');
    c=0;
    Xrep=0;
    n=size(X,2);
    while (any(sum(Xrep,1)==0)||any(sum(Xrep,2)==0))&&c<100
        if usingboot
            i=randi(n,1,500);
            Xrep=X(:,i);
        else
            Xrep=X(:,randperm(n));
            Xrep=Xrep(:,1:500);
        end
            c=c+1;
    end
	fprintf('......%d\n',c);
end
