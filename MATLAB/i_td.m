function [A0,A1]=i_td(XM0,XM1,methodid)

if nargin<3, methodid=1; end

switch methodid
    case 1
        Xhat0=do_td(XM0);
        Xhat1=do_td(XM1);
    case 2
        XM(:,:,:,1)=XM0;
        XM(:,:,:,2)=XM1;
        Xhat=do_td(single(XM));
        Xhat0=Xhat(:,:,:,1);
        Xhat1=Xhat(:,:,:,2);
end
A0=mean(Xhat0,3);
A1=mean(Xhat1,3);

%A0=A0./max(abs(A0(:)));
%A0=round(A0,5);
end

function Xhat0=do_td(XM0)
    T0=tensor(XM0);
    [~,U1]=cp_als(T0,5,'printitn',0);
    M2=cp_als(T0,5,'maxiters',100,...
              'init',U1,'printitn',0);
    %Use HOSVD initial guess
    %Use the 'nvecs' option to use the leading mode-n singular vectors as the initial guess.
    %M2 = cp_als(T0,5,'init','nvecs','printitn',0);
    fM0=full(M2);
    Xhat0=fM0.data;
end