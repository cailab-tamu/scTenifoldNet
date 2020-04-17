addpath('..\');
X0=readmatrix("X.txt");
X1=readmatrix("Y.txt");
t=readtable("genelist.txt","ReadVariableNames",false);
genelist=string(t.Var1);
tic;
T=sctenifoldnet_m(X0,X1,genelist,true);
toc;