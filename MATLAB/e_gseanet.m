function [G]=e_gseanet(Tf)
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
    nodenamesfull{k}=sprintf('%d_%s',k,Tf.pathway{k});
    a=sprintf('%d\\_%s',k,Tf.pathway{k});
    a=extractBefore(a,min(20,length(a)));
    nodenames{k}=a;
end

B=A.*(abs(A)>quantile(abs(A(:)),0.95));
% G=digraph(A,Tf.pathway);
G=digraph(B,nodenames);
LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
LWidths(LWidths==0)=1e-5;
%%
figure;
plot(G,'LineWidth',LWidths,'NodeLabel',...
    nodenames,'NodeFontAngle','normal',...
    'NodeFontSize',12);
p.MarkerSize = 7;
p.Marker = 's';
p.NodeColor = 'r';

%%
bins = conncomp(G);
for k=1:max(bins)
    fprintf('\nGroup %d ---------------\n',k);
    fprintf('%s\n',nodenamesfull{bins==k});    
end
fprintf('---------------\n');
