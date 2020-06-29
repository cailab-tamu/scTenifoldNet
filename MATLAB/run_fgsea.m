function [s]=run_fgsea(T)
if isempty(FindRpath)
   error('Rscript.ext is not found.');
end

oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/fgsea');
cd(pth);
if exist('output.txt','file'), delete('output.txt'); end
T.genelist=upper(T.genelist);
writetable(T,'input.txt');
RunRcode('script.R');
pause(1);
if exist('output.txt','file')
    s=readtable('output.txt',"Delimiter",',');
else
    s=[];
end
if exist('input.txt','file'), delete('input.txt'); end
if exist('output.txt','file'), delete('output.txt'); end
cd(oldpth);
