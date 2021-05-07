function [T]=i_dr(aln0,aln1,genelist,dosort)
% DR - differential regulatory gene identification 
    if nargin<4, dosort=true; end    
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
end



