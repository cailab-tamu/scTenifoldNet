function t = full(t)
%FULL Convert a symktensor to a symtensor.
%
%   T = FULL(C) converts a symktensor to a symtensor.
%
%   Examples
%   X = symktensor([3; 2], ones(4,2));
%   Y = full(A) %<-- equivalent dense tensor
%
%   See also SYMKTENSOR, TENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2015, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2015) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt

n = size(t,1);
m = ndims(t);
sz = nchoosek(m+n-1,m);

I = zeros(sz,m);
for loc = 1:sz
    if loc == 1
        I(loc,:) = ones(1,m);
    else
        I(loc,:) = I(loc-1,:);
        j = m;
        while (I(loc,j) == n)
            j = j - 1;
        end
        I(loc,j:m) = I(loc,j)+1;
    end
end

tnew = symtensor(@ones,m,n);
vals = entry(t,I);
tnew(I) = vals;
t = tnew;