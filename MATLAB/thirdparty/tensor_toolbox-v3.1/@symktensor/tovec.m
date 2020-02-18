function x = tovec(S,nolambda)
%TOVEC Convert symktensor to vector representation.
%
%   V = TOVEC(S) converts a symktensor to a vector. It stacks the LAMBDA
%   vector on top of a vectorized version of the matrix X.
%
%   V = TOVEC(S,TRUE) just returns a vectorized version of the matrix
%   X. It requires LAMBDA=1.
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

if exist('nolambda','var') && nolambda 
    if any(S.lambda ~= 1)
        error('Not all lambda values are 1.')
    end
    x = reshape(S.u,[],1);
else
    x = [S.lambda; reshape(S.u,[],1)];
end
