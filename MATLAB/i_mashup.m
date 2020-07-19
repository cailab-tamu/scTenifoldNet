function [aln0,aln1]=i_mashup(XM0,XM1,ndim)
    if nargin<3, ndim=30; end
    
    nsubsmpl=size(XM0,3);
    ngene=size(XM0,1);
    aln0=i_mashup_code(XM0,nsubsmpl,ngene,ndim)';
    aln1=i_mashup_code(XM1,nsubsmpl,ngene,ndim)';
end

function [x]=i_mashup_code(XM,nsubsmpl,ngene,ndim)
    RR_sum = zeros(ngene);
    for i = 1:nsubsmpl
      A = XM(:,:,i);
      if ~isequal(A,A'), A = A + A'; end
      % A = A + diag(sum(A, 2) == 0);      
      Q = rwr(A, 0.5);
      R = log(Q + 1/ngene); % smoothing
      RR_sum = RR_sum + R * R';
    end
    clear R Q A
    [V, d] = eigs(RR_sum, ndim);
    x = diag(sqrt(sqrt(diag(d)))) * V';
end
