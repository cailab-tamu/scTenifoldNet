function sz = size(t,idx)
%SIZE Size of symktensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = SIZE(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   See also SYMKTENSOR, SYMKTENSOR/NDIMS.
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


if isempty(t.lambda)
    sz = [];
end

if exist('idx','var')
    sz = size(t.u, 1);
else
    sz = size(t.u, 1) * ones(1,t.m);  
end
