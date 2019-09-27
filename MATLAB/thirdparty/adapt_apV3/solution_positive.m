%load GeneFindingProblem.mat;
% offset=0; % Change this prefence offset to change sensitivity
% [idx,netsim,dpsim,expref]=apcluster(s,p+offset,'maxits',15,'sparse');
% ap_exon=(idx(1:end-1)~=75067); % remove non-exon exemplar and identify exons
% true_positive_rate=sum(ap_exon.*refseq_exon)/sum(refseq_exon)
% false_positive_rate=sum(ap_exon.*refseq_intron)/sum(refseq_intron)

function [true_positive,false_positive] = solution_positive(refseq_exon,refseq_intron,labelid,classend,Sid)
% remove non-exon exemplar and identify exons
ap_exon = (labelid(1:end-1,Sid)~=classend); 
true_positive=sum(ap_exon.*refseq_exon)/sum(refseq_exon);
false_positive=sum(ap_exon.*refseq_intron)/sum(refseq_intron);