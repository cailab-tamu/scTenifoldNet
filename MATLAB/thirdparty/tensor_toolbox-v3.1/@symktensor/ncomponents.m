function n = ncomponents(t)
%NCOMPONENTS Number of components for a symktensor.
%
%   NCOMPONENTS(T) returns the number of components in the symktensor T.
%   This is the size of the lambda vector, equivalently the number of
%   columns in the factor matrix.
%
%   S = symktensor(3, symtensor(@rand,4,3)); %<--Random symktensor
%   ncomponents(S) %<--Returns 3
%
%   See also SYMKTENSOR
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


n = length(t.lambda);
