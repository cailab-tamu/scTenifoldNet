function res = innerprod(X,Y)
%INNERPROD Efficient inner product with a sumtensor.
%
%   R = INNERPROD(X,Y) efficiently computes the inner product between
%   two tensors X and Y, where X is a sumtensor.
%
%   See also TENSOR/INNERPROD, SPTENSOR/INNERPROD, TTENSOR/INNERPROD, 
%   KTENSOR/INNERPROD
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

% X is a sumtensor

res = innerprod(X.part{1},Y);
for i = 2:length(X.part)
    res = res + innerprod(X.part{i},Y); 
end



