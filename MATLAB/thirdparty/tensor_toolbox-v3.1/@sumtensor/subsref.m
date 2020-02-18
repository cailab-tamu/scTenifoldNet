function B = subsref(A, S) 
%SUBSREF Subscript reference for sumtensor.
%
%   T.part{i} returns the ith part of the sumtensor T
%
%   Examples
%   T1 = tensor(rand(3,3,3));
%   T2 = sptensor([1 1 1; 3 1 2; 1 1 3], 1, [3,3,3]);
%   T = sumtensor(T1,T2); 
%   T.part{2} %<--Returns the symmetric tensor
%
%   See also SUMTENSOR.
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

B = builtin('subsref', A, S);

