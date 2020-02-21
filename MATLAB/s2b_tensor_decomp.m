% load examplenetworks.mat A0 A1
% clear X
for k=1:length(Ax0)
    Xx0(:,:,k)=full(Ax0{k});
end
for k=1:length(Ax1)
    Xx1(:,:,k)=full(Ax1{k});
end

%%
R=5;
[A,B,C,lam]=cp3_DCPD(Xx0,R);
A0=tensor_reconst(A,B,C,lam);
[A,B,C,lam]=cp3_DCPD(Xx1,R);
A1=tensor_reconst(A,B,C,lam);


A0=mean(A0,3);
A0=A0./max(abs(A0(:)));
A0=round(A0,1);

A1=mean(A1,3);
A1=A1./max(abs(A1(:)));
A1=round(A1,1);

%%

function [A,B,C,lam]=cp3_DCPD(X,R)
Iter1=50; Iter2=3;
A=randn(size(X,1),R);B=randn(size(X,2),R);C=randn(size(X,3),R);
for i=1:Iter1
    for r=1:R
        ind=setdiff(1:R,r);
        lam=lam_find(X,A(:,ind),B(:,ind),C(:,ind));
        EST=tensor_reconst(A(:,ind),B(:,ind),C(:,ind),lam);
        for i2=1:Iter2
        A(:,r)=tens2mat(X-EST,1)/(kr(B(:,r),C(:,r))');
        B(:,r)=tens2mat(X-EST,2)/(kr(A(:,r),C(:,r))');
        C(:,r)=tens2mat(X-EST,3)/(kr(A(:,r),B(:,r))');
        end 
    end
end
lam=lam_find(X,A,B,C);
end

%%%%%% FUNCTIONS %%%%%%%%%
function Xhat=tensor_reconst(A,B,C,lam)
    Xhat=(A*diag(lam))*(kr(C,B)');
    Xhat=reshape(Xhat,[size(A,1),size(B,1),size(C,1)]);
end
 %%%%%%%%%%%%%%%%%%%%%%
function lam=lam_find(X,A,B,C)
    for r=1:size(A,2)
       D(:,r)=vec(tensor_reconst(A(:,r),B(:,r),C(:,r),1));
    end
    lam=pinv(D)*X(:);
end


 %%%%%%%%%%%%%%%%%%%%%%
function [X_mat]=tens2mat(X,mode)
    ORDERS=[1 3 2;2 3 1;3 2 1];
    X_mat=reshape(permute(X,ORDERS(mode,:)),size(X,mode),numel(X)/size(X,mode));
end

 %%%%%%%%%%%%%%%%%%%%%%
function v = vec(x)
    v = reshape(x,numel(x),1);
end
 %%%%%%%%%%%%%%%%%%%%%%
function AB = kr(A,B)
[I,F]=size(A);[J,~]=size(B);AB=zeros(I*J,F);
for f=1:F
    AB(:,f)=vec(B(:,f)*A(:,f).');
end
end

