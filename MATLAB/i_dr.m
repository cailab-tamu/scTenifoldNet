function [T]=i_dr(aln0,aln1,genelist,doplot,dosort)
% DR - differential regulatory gene identification 
    if nargin<5, dosort=true; end
    if nargin<4, doplot=false; end
    if nargin<3, genelist=string(num2cell(1:size(aln0,1)))'; end
    drdist=vecnorm(aln0-aln1,2,2).^2;
    drdist=drdist./norm(drdist);
    FC=drdist./mean(drdist);
    pValues=chi2cdf(FC,1,'upper');
    pAdjusted = mafdr(pValues,'BHFDR',true);
    if size(genelist,1)==1, genelist=genelist'; end
    sortid=(1:length(genelist))';
    if size(genelist,2)>1, genelist=genelist'; end
    T=table(sortid,genelist,drdist,FC,pValues,pAdjusted);
    if dosort
        T = sortrows(T,'drdist','descend');
    end
    if doplot
        pd = makedist('Gamma','a',0.5,'b',2);
        qqplot(FC,pd);
        [~,i]=sort(FC);
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist(i)};
    end
end

function txt = i_myupdatefcn1(~,event_obj,g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
end

