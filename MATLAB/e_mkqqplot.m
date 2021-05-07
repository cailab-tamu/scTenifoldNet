function e_mkqqplot(T)

pd = makedist('Gamma','a',0.5,'b',2);
qqplot(T.FC,pd);
[~,i]=sort(T.FC);
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1,T.genelist(i)};
    
end

function txt = i_myupdatefcn1(~,event_obj,g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
end

