function [s]=e_fgsea(T,rmribo)
if nargin<2
    rmribo=false;
end
if isempty(FindRpath)
   error('Rscript.ext is not found.');
end

oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/fgsea');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

if exist('output.txt','file'), delete('output.txt'); end
T.genelist=upper(string(T.genelist));
if rmribo
    [gribo]=get_ribosomalgenes;
    i=~ismember(T.genelist,gribo);
    T=T(i,:);
end
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
