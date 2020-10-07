function m = size(t,idx)
%SIZE Size of a sumtensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = size(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   See also SUMTENSOR, SUMTENSOR/NDIMS.
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

if isempty(t.part)
    m = [];
elseif exist('idx','var')
    m = size(t.part{1},idx);
else
    m = size(t.part{1});
end
    
