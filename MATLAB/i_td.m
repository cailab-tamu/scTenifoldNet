function A0=i_td(XM0,methodid)

if nargin<1
    methodid=2;
end

switch methodid
    case 1
    R=5;
    [A,B,C,lam]=cp3_DCPD(XM0,R);
    Xhat0=tensor_reconst(A,B,C,lam);
    
    case 2
    T0=tensor(XM0);
    [~,U1]=cp_als(T0,5,'printitn',0);
    M2=cp_als(T0,5,'maxiters',100,'init',U1,'printitn',0);
    fM0=full(M2);
    Xhat0=fM0.data;
end


A0=mean(Xhat0,3);
%A0=A0./max(abs(A0(:)));
%A0=round(A0,5);
