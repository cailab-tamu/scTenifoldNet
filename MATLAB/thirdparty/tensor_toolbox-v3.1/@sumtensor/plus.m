function t = plus( t, x )
%PLUS Plus for sumtensor.
%
%   T = T + A adds A to the sumtensor parts, where T is a sumtensor and 
%   A is any valid sumtensor component.
%
%   T = T + {A,B,C} adds A, B, and C to the sumtensor parts, assuming they
%   are valid sumtensor components.
%
%   Note that new parts are appended to the end of T, even if A + T is called
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
% The full license terms can be found in the file LICENSE.t


    
if iscell(x)
    t = sumtensor(t.part{:}, x{:});
else    
    t = sumtensor(t.part{:}, x);
end
