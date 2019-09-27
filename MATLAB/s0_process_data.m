% [X,genelist,celllist]=sc_readtsvfile("expression_Aging_mouse_brain_portal_data_updated.txt");
% celllist(1)
% celllist(1)=[];
% X=sparse(X);
% save expression_Aging_mouse_brain_portal_data_updated X genelist celllist

% T=readtable('meta_Aging_mouse_brain_portal_data.txt','delimiter','\t');



sum(string(T.cell_type)=="MG")
sum(string(T.cell_type_by_age)=="MG_2-3mo")
sum(string(T.cell_type_by_age)=="MG_21-22mo")
i0=string(T.cell_type_by_age)=="MG_2-3mo";
i1=string(T.cell_type_by_age)=="MG_21-22mo";
X0=X(:,i0);
X1=X(:,i1);

i=string(T.cell_type)=="MG";
agegroup_mg=T.all_cells_by_age(i);
Xmg=full(X(:,i));


%%

sum(string(T.cell_type_by_age)=="NEUR_mature_2-3mo")
sum(string(T.cell_type_by_age)=="NEUR_mature_21-22mo")

i=string(T.cell_type)=="NEUR_mature";
agegroup_ne=T.all_cells_by_age(i);
Xne=full(X(:,i));

s1=run_phate(Xne,3);
figure; i_myscatter(s1,grp2idx(agegroup));
s2=sc_tsne(Xne,3,false,false);
figure; i_myscatter(s2,grp2idx(agegroup));

