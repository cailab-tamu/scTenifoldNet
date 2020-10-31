function [A0,A1]=i_td(XM0,XM1,methodid)

if nargin<3, methodid=1; end
if exist('tensor.m','file')~=2
    pw0=pwd;
    pw1=fileparts(which(mfilename));
    cd(pw1);
    addpath('thirdparty/tensor_toolbox');
    cd(pw0);
end
if exist('tensor.m','file')~=2
    error('Need thirdparty/tensor_toolbox');
end

switch methodid
    case 1
        Xhat0=do_td_cp(XM0);
        Xhat1=do_td_cp(XM1);
    case 2
        Xhat0=do_td_tucker(XM0);
        Xhat1=do_td_tucker(XM1);
    case 3
        XM(:,:,:,1)=XM0;
        XM(:,:,:,2)=XM1;
        Xhat=do_td_cp(single(XM));
        Xhat0=Xhat(:,:,:,1);
        Xhat1=Xhat(:,:,:,2);
end
A0=mean(Xhat0,3);
A1=mean(Xhat1,3);

A0=e_transf(A0);
A1=e_transf(A1);
end

function Xhat0=do_td_cp(XM0)
    T0=tensor(XM0);
    % Use HOSVD initial guess    
    M = cp_als(T0,5,'init','nvecs','printitn',1,'maxiters',200);
    % Xhat0=double(M);
    fM0=full(M);
    Xhat0=fM0.data;    
end

function Xhat0=do_td_tucker(XM0)
    T0=tensor(XM0);
    M=tucker_als(T0,5,'printitn',1);
    % Xhat0=double(M);
    fM0=full(M);
    Xhat0=fM0.data;
end
