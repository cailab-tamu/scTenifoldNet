function C = mtimes(A,B)
%MTIMES tensor-scalar multiplication.
% 
%   C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   symtensor and the other argument is a scalar.
% 
%   For symtensor-symtensor array multiplication, use TIMES or 'A .* B'.
% 
%   Examples
%   X = symtenrand([4,4,4])
%   W = 5 * X
%
%   See also SYMTENSOR, SYMTENSOR/TIMES
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


%%
if isscalar(B)
    C = A;
    C.val = B * C.val;
    return;
end

if isscalar(A)
    C = B;
    C.val = A * C.val;
    return;
end

error('Mtimes only supports a symtensor times a scalar');





