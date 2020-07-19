function Q = rwr(A, p)
if nargin<2, p=0.5; end
  n = size(A, 1);
  A = A - diag(diag(A));
  % A = A + diag(sum(A) == 0);  
  B = A./sum(A);  
  Q=mldivide(eye(n)-(1-p)*B, p*eye(n));
end
