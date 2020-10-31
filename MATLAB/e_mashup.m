function [aln0,aln1]=e_mashup(XM0,XM1,ndim)
% Mashup for obtaining high-quality, compact topological 
% feature representations of genes from one or more interaction
% networks constructed from heterogeneous data types.
%
% ref: doi: 10.1016/j.cels.2016.10.017
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
      if ~issymmetric(A), A = A + A'; end
      % if ~isequal(A,A'), A = A + A'; end
      % % A = A + diag(sum(A, 2) == 0);
      Q = rwr(A, 0.5);
      R = log(Q + 1/ngene); % smoothing
      RR_sum = RR_sum + R * R';
    end
    clear R Q A
    [V, d] = eigs(RR_sum, ndim);
    x = diag(sqrt(sqrt(diag(d)))) * V';
end

function Q = rwr(A, p)
if nargin<2, p=0.5; end
  n = size(A, 1);
  A = A - diag(diag(A));
  % A = A + diag(sum(A) == 0);  
  B = A./sum(A);  
  Q=mldivide(eye(n)-(1-p)*B, p*eye(n));
end
