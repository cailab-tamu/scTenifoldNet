function sz = size(t,idx)
%SIZE Dimensions of a symmetric tensor.
%  
%   D = SIZE(T) returns the size of a symtensor T. 
%
%   I = SIZE(T,DIM) returns the size of symtensor T in the dimension 
%   specified by the scalar DIM.
%
%   See also SYMTENSOR, SYMTENSOR/NDIMS, SIZE.
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


if (t.m == 0)
    sz = [];
elseif exist('idx','var')
    if 1<=idx && idx<=t.m  %Bounds check
        sz = t.n;
    else
        error('Index exceeds tensor dimensions');
    end
else
    sz = t.n * ones(1,t.m);  
end
