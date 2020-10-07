function n = ncomponents(t)
%NCOMPONENTS Number of components for a ktensor.
%
%   NCOMPONENTS(T) returns the number of components in the ktensor T.
%
%   X = ktensor(ones(4,1), rand(2,4), randn(3,4), randi(5,4,4));
%   ncomponents(X) %<--Returns 4
%
%   See also KTENSOR
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
