function [Tp,Tn]=run_fgsea2(genelist,generank,dbfile)
oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/fgsea');
cd(pth);

if isempty(FindRpath)
   cd(oldpth);
   error('Rscript.exe is not found.');
end

if exist('output.txt','file'), delete('output.txt'); end
if nargin<2 || isempty(generank), generank=1:length(genelist); end
if nargin<3, dbfile='bp'; end

if all(generank== floor(generank))
    generank=generank./1.1;
end

fid=fopen('input.txt','w');
fprintf(fid,'ID\tt\n');
% https://github.com/ctlab/fgsea/blob/master/inst/extdata/naive.vs.th1.rnk
for k=1:length(genelist)
    fprintf(fid,'%s\t%f\n',upper(genelist(k)),generank(k));
end
fclose(fid);
switch lower(dbfile)
    case 'bp'
        RunRcode('script_bp.R');
    case 'mf'
        RunRcode('script_mf.R');
end


T=readtable('output.txt');
% T=sortrows(T,3);
T=sortrows(T,'pval','ascend');
Tn=T(T.ES<0,:);
Tp=T(T.ES>0,:);
if exist('input.txt','file'), delete('input.txt'); end
if exist('output.txt','file'), delete('output.txt'); end
readtable('output_top20.txt')
if exist('output_top20.txt','file'), delete('output_top20.txt'); end
cd(oldpth);
