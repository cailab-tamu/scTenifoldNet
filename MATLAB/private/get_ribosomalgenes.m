function [g]=get_ribosomalgenes
options = weboptions('Timeout',21);
websave('a.txt','https://www.genenames.org/cgi-bin/genegroup/download?id=1054&type=branch',options);
% t=readtable('a.txt','PreserveVariableNames',true);
warning off 
t=readtable('a.txt','Range','A:B');
warning on
g=string(t.ApprovedSymbol);
delete('a.txt');
