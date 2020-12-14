function e_fgseanet(Tf,jaccd)
% Merge similar gene sets (Jaccard index > 0.8) in fGSEA report
if nargin<2, jaccd=0.8; end
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
B=A.*(A>jaccd);
% G=digraph(A,Tf.pathway);
G=digraph(B,nodenames);
LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
LWidths(LWidths==0)=1e-5;
%%
figure;
p=plot(G,'NodeLabel',nodenames,'NodeLabelMode','auto'); 

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

tmpName=[tempname,'.txt'];
fid=fopen(tmpName,'w');
for k=1:max(bins)
    fprintf(fid,'\nEnriched Function Group %d\n',k);
    fprintf(fid,'\t%s\n',nodenamesfull{bins==idx(k)});    
end
%fprintf(fid,'---------------\n');
fclose(fid);
[status]=system(['notepad "' tmpName '" &']);
if status~=0
   edit(tmpName);
end

