function D = L2_distance(a, b, df)
%
% D = L2_distance(a, b, df=0);
%
% a and b are k*n matrices, representing n k-dimensional points
% D(i,j) is the L2 distance between a(i) and b(j)

  if nargin == 2, df = 0; end

  if (size(a,1) == 1)
    a = [a; zeros(1,size(a,2))]; 
    b = [b; zeros(1,size(b,2))]; 
  end
  
  aa = sum(a.*a);
  bb = sum(b.*b);
  ab = a'*b; 
  D = real(sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
  
  if (df==1)
    D = D.*(1-eye(size(D)));
  end
