function Z = plus( X, Y )
%PLUS Binary addition (+) for symtensors.
%
%   PLUS(A,B) is called for the syntax 'X + Y' when X and Y are symtensors.
%   A and B must be the same size, unless one is a scalar. A scalar can be
%   added to a symtensor of any size.
%
%   See also SYMTENSOR, TENSOR/PLUS.
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

Z = tenfun(@plus,X,Y);


    


