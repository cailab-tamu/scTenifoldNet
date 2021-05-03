function [OUT]=e_fgseanet(Tf,varargin)
% Merge similar gene sets (Jaccard index > cutoff) in fGSEA report

   p = inputParser;
   addOptional(p,'PlotNetwork',false,@islogical);
   addOptional(p,'JaccardCutoff',0.5,@(x) x>0 & x<1)
   parse(p,varargin{:});
   plotnetwork=p.Results.PlotNetwork;
   jaccardcutoff=p.Results.JaccardCutoff;

if size(Tf,1)<5
    error('fGSEA output table is too short.')
end

n=size(Tf.leadingEdge,1);
A=zeros(n);
for i=1:n-1
    for j=i+1:n
        a=strsplit(Tf.leadingEdge{i},";");
        b=strsplit(Tf.leadingEdge{j},";");        
        A(i,j)=length(intersect(a,b))./length(unique(union(a,b)));
        A(j,i)=A(i,j);
    end
end
%%
nodenames=Tf.pathway;
nodenamesfull=Tf.pathway;
for k=1:n
    % nodenamesfull{k}=sprintf('%d_%s',k,Tf.pathway{k});
    % nodenamesfull{k}=sprintf('%s',Tf.pathway{k});
    nodenamesfull{k}=sprintf('%d_%s',k,Tf.pathway{k});
    %a=sprintf('%d\\_%s',k,Tf.pathway{k});
    %a=extractBefore(a,min(20,length(a)));
    nodenames{k}=sprintf('%d',k);
end
%%
%B=A.*(abs(A)>quantile(abs(A(:)),0.95));
B=A.*(A>jaccardcutoff);
% G=digraph(A,Tf.pathway);
G=digraph(B,nodenames);
LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
LWidths(LWidths==0)=1e-5;

if plotnetwork
figure;
p=plot(G,'NodeLabel',nodenames,'NodeLabelMode','auto'); 
end

% p=plot(G,'NodeLabel',nodenames,'NodeFontAngle','normal',...
%     'NodeFontSize',12);
% if ~isempty(LWidths)
%     p.LineWidth=LWidths;
% end
% p.MarkerSize = 7;
% p.Marker = 's';
% p.NodeColor = 'r';

%%
[bins,binsizes] = conncomp(G);
[~,idx]=sort(binsizes,'descend');

OUT=cell(max(bins),2);

tmpName=[tempname,'.txt'];
fid=fopen(tmpName,'w');
for k=1:max(bins)
    fprintf(fid,'\nEnriched Function Group %d\n',k);
    vi=find(bins==idx(k));
    Gx=[];
    for kk=1:length(vi)
        a=strsplit(Tf.leadingEdge{vi(kk)},";");
        Gx=[Gx a];
    end
    Gx=unique(Gx,'stable');
    OUT{k,1}=string(Gx);
    
    fprintf(fid,'%s ',string(Gx));
    fprintf(fid,'\n');
    fprintf(fid,'\t%s\n',nodenamesfull{bins==idx(k)});
    OUT{k,2}=deblank(sprintf('%s\n',nodenamesfull{bins==idx(k)}));
end

%fprintf(fid,'---------------\n');
fclose(fid);
[status]=system(['notepad "' tmpName '" &']);
if status~=0
   edit(tmpName);
end
end
