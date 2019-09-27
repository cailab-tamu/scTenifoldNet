if algorithm == 1
  [Smax, Sid] = max(Sil);
  NCopt = NC(Sid);
  fprintf('\n## Clustering solution by adaptive Affinity Propagation:\n');
  fprintf('  Optimal number of clusters is %d, Silhouette = %g,\n',NCopt,Smax);
  if Smax < 0.3
      R = length(NC);
      R = ceil(R/2):R;
      [Tmax, Q] = max(Silmin(R));
      Sid = R(Q);
      NCopt2 = NC(Sid);
      fprintf('  If Silhouette values are small & NCs are large, Optimal NC is %d,\n',NCopt2);
      fprintf('  where min Silhouette of single cluster is %g.\n',Tmax);
      fprintf('  The optimal solution (class labels) is in labels(:,Sid)');
  end
  fprintf('\n## Silhouette values at different NCs: [NC;Sil;Silmin] \n');
  disp([NC;Sil;Silmin]);
  if id == 14
    [TP,FP] = solution_positive(refseq_exon,refseq_intron,labelid,nsubset,Sid);
    fprintf('\n## exon identification: true positive rate %g, false positive rate %g\n',TP,FP);
  end
 
elseif  algorithm == 2
  [Smax, Sid] = max(Sil);
  NCopt = NC(Sid);
  fprintf('\n## Clustering solution searched by Affinity Propagation:\n');
  fprintf('  Optimal number of clusters is %d, Silhouette = %g,\n',NCopt,Smax);
  if Smax < 0.3
      R = length(NC);
      R = ceil(R/2):R;
      [Tmax, Q] = max(Silmin(R));
      Sid = R(Q);
      NCopt2 = NC(Sid);
      fprintf('  If Silhouette values are small & NCs are large, Optimal NC is %d,\n',NCopt2);
      fprintf('  where min Silhouette of single cluster is %g.\n',Tmax);
      fprintf('  The optimal solution (class labels) is in labels(:,Sid)');
  end
  fprintf('\n## Silhouette values at different NCs: [NC;Sil;Silmin] \n');
  disp([NC;Sil;Silmin]);
  if id == 14
    [TP,FP] = solution_positive(refseq_exon,refseq_intron,labelid,nsubset,Sid);
    fprintf('\n## exon identification: true positive rate %g, false positive rate %g\n',TP,FP);
  end
    
else
    NCs = unique(labels);
    NCopt = length(NCs);
    Sid = 1;
    [C, labels] = ind2cluster(labels);
    [NC,Sil,Silmin] = solution_evaluation(data,M,labels,NCopt,...
      iend,simatrix,nrow,type,cut);
    fprintf('$$ Clustering solution by original Affinity Propagation:\n');
    fprintf('  Optimal number of clusters is %d, Silhouette = %g,\n',NCopt,Sil);
    fprintf('  where min Silhouette of single cluster is %g.\n',Silmin);
    fprintf('  The optimal solution (class labels) is in labels(:,Sid)');
    if id == 14
      [TP,FP] = solution_positive(refseq_exon,refseq_intron,labels,labels(end),Sid);
      fprintf('\n## exon identification: true positive rate %g, false positive rate %g\n',TP,FP);
    end
end