function  [XM0,XM1]=i_nc(X0,X1,nsubsmpl,ncom,usebootstrp,csubsmpl)

if nargin<6, csubsmpl=500; end       % number of cells in subsamples
if nargin<5, usebootstrp=false; end   % using m-out-of-n bootstrap (false by default)
                                    % using jackknife (by default)
if nargin<4, ncom=3; end    % number of components
if nargin<3, nsubsmpl=10; end      % number of subsamples 

    n=size(X0,1);
    XM0=zeros(n,n,nsubsmpl);
    XM1=zeros(n,n,nsubsmpl);
    for k=1:nsubsmpl
        fprintf('network...%d of %d\n',k,nsubsmpl);
        
        n0=size(X0,2);
        if usebootstrp % bootstrap 
            i=randi(n0,1,csubsmpl);
            Xrep=X0(:,i);
        else     % jackknife
            Xrep=X0(:,randperm(n0));
            Xrep=Xrep(:,1:csubsmpl);
        end    
        
        %if any(sum(Xrep,1)==0)||any(sum(Xrep,2)==0)
        %    Xrep=i_goodrep(X0,usingboot);
        %end
        
        A=sc_pcnetpar(Xrep,ncom,false);
        A=A./max(abs(A(:)));        
        XM0(:,:,k)=A.*(abs(A)>quantile(abs(A(:)),0.95));

        n1=size(X1,2);        
        if usebootstrp % bootstrap 
            i=randi(n1,1,csubsmpl);
            Xrep=X1(:,i);
        else     % jackknife
            Xrep=X1(:,randperm(n1));
            Xrep=Xrep(:,1:csubsmpl);
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
            i=randi(n,1,csubsmpl);
            Xrep=X(:,i);
        else
            Xrep=X(:,randperm(n));
            Xrep=Xrep(:,1:csubsmpl);
        end
            c=c+1;
    end
	fprintf('......%d\n',c);
end
