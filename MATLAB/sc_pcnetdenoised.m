function [A]=sc_pcnetdenoised(X,varargin)
%Construct network using scTenifoldNet (accurate, robust, but slow)
% A=sc_pcnetdenoised(X);
%
% X is gene x cell matrix
% 

if exist('sctenifoldnet','file')~=2
    error('Requires sctenifoldnet.m');
end

    if nargin<1
        error(sprintf('USAGE: A=sc_pcnetdenoised(X);\n       A=sc_pcnetdenoised(X,''smplmethod'',''Jackknife'');'));
    end
    
   p = inputParser;
   addOptional(p,'smplmethod',"Jackknife",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["jackknife","bootstrap"]));
   addOptional(p,'tdmethod',"CP",@(x) (isstring(x)|ischar(x))&ismember(upper(string(x)),["CP","TUCKER"]));
   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
   addOptional(p,'savegrn',false,@islogical);
   parse(p,varargin{:});
   tdmethod=p.Results.tdmethod;
   nsubsmpl=p.Results.nsubsmpl;
   csubsmpl=p.Results.csubsmpl;
   smplmethod=p.Results.smplmethod;
   savegrn=p.Results.savegrn;
   
   switch upper(tdmethod)
       case "CP"
           tdmethod=1;
       case "TUCKER"
           tdmethod=2;
   end
   switch lower(smplmethod)
       case "jackknife"
           usebootstrp=false;
       case "bootstrap"
           usebootstrp=true;
   end
    pw0=pwd;    
    pw1=fileparts(mfilename('fullpath'));
    cd(pw1);
    pth=fullfile(pw1,'thirdparty','tensor_toolbox');
    addpath(pth);
    cd(pw0);
    
    if exist('tensor.m','file')~=2
        error('Need thirdparty/tensor_toolbox');
    end
    if exist('sc_pcnetpar.m','file')~=2
        error('Need sc_pcnetpar.m in scGEAToolbox https://github.com/jamesjcai/scGEAToolbox');
    end
    
    
    X=sc_norm(X,"type","libsize");
    %X=sc_transform(X);
    %tic
    [XM]=i_nc(X,nsubsmpl,3,csubsmpl,usebootstrp);
    %toc

    %tic
    disp('Tensor decomposition')
    [A]=i_td1(XM,tdmethod);
    %toc
    if savegrn
        tstr=matlab.lang.makeValidName(datestr(datetime));
        save(sprintf('A_%s',tstr),'A','-v7.3');
    end    
end