function T=sctenifoldnet_m(X0,X1,genelist,doplot)
    if nargin<3
        error('USAGE: T=sctenifoldnet_m(X0,X1,genelist)');
    end
    if nargin<4, doplot=false; end
    if size(X0,1)~=size(X1,1)
        error('X0 and X1 need the same number of rows.');
    end
    if size(X0,1)~=length(genelist)
        error('Length of genelist should be the same as the number of rows of X0 or X1.');
    end
    pw0=pwd;
    pw1=fileparts(which(mfilename));
    cd(pw1);
    addpath('thirdparty\tensor_toolbox-v3.1\');
    cd(pw0);
    if exist('sc_pcnet.m','file')~=2
        error('scGEAToolbox installation is required. Link: https://github.com/jamesjcai/scGEAToolbox');
    end    
    
    X0=sc_norm(X0,"type","libsize");
    X1=sc_norm(X1,"type","libsize");
    
    [XM0,XM1]=i_nc(X0,X1);
    [A0,A1]=i_td(XM0,XM1);
    [aln0,aln1]=i_ma(A0,A1);
    T=i_dr(aln0,aln1,genelist,doplot);
end

