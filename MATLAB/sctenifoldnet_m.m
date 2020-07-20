function T=sctenifoldnet_m(X0,X1,genelist,varargin)
    % X0 and X1 are gene x cell matrices
    if nargin<2
        error(sprintf('USAGE: T=sctenifoldnet_m(X0,X1);\n       T=sctenifoldnet_m(X0,X1,genelist,''qqplot'',true);'));
    end
    if nargin<3, genelist=string(num2cell(1:size(X0,1)))'; end
    
   p = inputParser;
   addOptional(p,'qqplot',false,@islogical);
   addOptional(p,'smplmethod',"Jackknife",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["jackknife","bootstrap"]));
   addOptional(p,'tdmethod',"CP",@(x) (isstring(x)|ischar(x))&ismember(upper(string(x)),["CP","TUCKER"]));
   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
   parse(p,varargin{:});
   doqqplot=p.Results.qqplot;
   tdmethod=p.Results.tdmethod;
   nsubsmpl=p.Results.nsubsmpl;
   smplmethod=p.Results.smplmethod;
   
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
    pw1=fileparts(which(mfilename));
    cd(pw1);
    addpath('thirdparty\tensor_toolbox-v3.1\');
    cd(pw0);
    if exist('sc_pcnet.m','file')~=2
        error('Need sc_pcnet.m in the scGEAToolbox https://github.com/jamesjcai/scGEAToolbox');
    end    
    
    X0=sc_norm(X0,"type","libsize");
    X1=sc_norm(X1,"type","libsize");
    
    [XM0,XM1]=i_nc(X0,X1,nsubsmpl,3,usebootstrp);
    
    [A0,A1]=i_td(XM0,XM1,tdmethod);
    A0=0.5*(A0+A0');
    A1=0.5*(A1+A1');
    [aln0,aln1]=i_ma(A0,A1);
    
    % [aln0,aln1]=i_mashup(XM0,XM1);
    T=i_dr(aln0,aln1,genelist,doqqplot);
end

