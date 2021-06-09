function [T,A0,A1]=sctenifoldnet(X0,X1,genelist,varargin)
% T=sctenifoldnet_m(X0,X1,genelist);
%
% X0 and X1 are gene x cell matrices
% 
    if nargin<2
        error(sprintf('USAGE: T=sctenifoldnet(X0,X1);\n       T=sctenifoldnet(X0,X1,genelist,''qqplot'',true);'));
    end
    if nargin<3, genelist=string(num2cell(1:size(X0,1)))'; end
    
   p = inputParser;
   addOptional(p,'qqplot',false,@islogical);
   addOptional(p,'smplmethod',"Jackknife",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["jackknife","bootstrap"]));
   addOptional(p,'tdmethod',"CP",@(x) (isstring(x)|ischar(x))&ismember(upper(string(x)),["CP","TUCKER"]));
   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
   addOptional(p,'savegrn',true,@islogical);
   parse(p,varargin{:});
   doqqplot=p.Results.qqplot;
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
   
    if size(X0,1)~=size(X1,1)
        error('X0 and X1 need the same number of rows.');
    end
    if size(X0,1)~=length(genelist)
        error('Length of genelist should be the same as the number of rows of X0 or X1.');
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
    if exist('sc_pcnet.m','file')~=2
        error('Need sc_pcnet.m in scGEAToolbox https://github.com/jamesjcai/scGEAToolbox');
    end    
    
    
    X0=sc_norm(X0,"type","libsize");
    X1=sc_norm(X1,"type","libsize");
    %X0=sc_transform(X0);
    %X1=sc_transform(X1);
    
    rng('default');
    
    tic
    disp('Sample 1/2 ...')
    [XM]=i_nc(X0,nsubsmpl,3,csubsmpl,usebootstrp);
    toc
    tic
    disp('Tensor decomposition')
    [A0]=i_td1(XM,tdmethod);
    [A0]=e_filtadjc(A0);
    toc
    if savegrn
        tic
        disp('Saving scGRN network 1/2')
        tstr=matlab.lang.makeValidName(datestr(datetime));
        save(sprintf('A0_%s',tstr),'A0','genelist','-v7.3');
        toc
    end
    tic
    disp('Sample 2/2 ...')
    [XM]=i_nc(X1,nsubsmpl,3,csubsmpl,usebootstrp);
    toc
    tic
    disp('Tensor decomposition')
    A1=i_td1(XM,tdmethod);
    A1=e_filtadjc(A1);
    toc
    if savegrn
        tic
        disp('Saving scGRN network 2/2')
        save(sprintf('A1_%s',tstr),'A1','genelist','-v7.3');
        toc
    end
    
    A0sym=0.5*(A0+A0');
    A1sym=0.5*(A1+A1');
    tic
    disp('Manifold alignment')
    [aln0,aln1]=i_ma(A0sym,A1sym);
    toc
    
    T=i_dr(aln0,aln1,genelist);
    
    if doqqplot
        e_mkqqplot(T);
    end
end
