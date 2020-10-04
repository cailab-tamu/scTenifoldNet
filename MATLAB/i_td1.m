function [A0]=i_td1(XM0,methodid)

if nargin<2, methodid=1; end
    
pw0=pwd;
pw1=fileparts(which(mfilename));
cd(pw1);
addpath('thirdparty/tensor_toolbox');
cd(pw0);
    
switch methodid
    case 1
        Xhat0=do_td_cp(XM0);        
    case 2
        Xhat0=do_td_tucker(XM0);    
end
A0=mean(Xhat0,3);
A0=e_transf(A0);
end




