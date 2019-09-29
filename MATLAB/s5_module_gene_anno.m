
gg=[genelist;genelist];
for k=1:max(C)
    i=C==k;
    if sum(i)>5
        fprintf('%d %d ',k,sum(i));
        fprintf('%s ',gg(i))
        fprintf('\n');
    end
end
%%
P_v=[];
G0_v=[]; G1_v=[];
J_v=[];


n=size(A0,1);

for k=1:max(C)
    
idx=find(C==k);
if length(idx)>20
    gid=idx;
    g0id=gid(gid<=n);
    g1id=gid(gid>n)-n;
    jacx=length(intersect(g0id,g1id))/length(union(g0id,g1id));
    if length(g0id)>length(g1id)
        J_v=[J_v; jacx];
    else
        J_v=[J_v; -jacx];
    end
    if ~isempty(g0id) && ~isempty(g1id)
        B0 = reshape(A0(g0id,g0id),[],1);
        B1 = reshape(A1(g1id,g1id),[],1);    
        [~,p]=kstest2(full(B0),full(B1));
        P_v=[P_v;p];
    else
        P_v=[P_v;nan];
    end
    G0_v=[G0_v;string(sprintf('%s ',genelist(g0id)))];
    G1_v=[G1_v;string(sprintf('%s ',genelist(g1id)))];
end
end

absJ_v=abs(J_v);
Tres=table(absJ_v,J_v,P_v,G0_v,G1_v);
writetable(Tres,'res_v','filetype','spreadsheet');
