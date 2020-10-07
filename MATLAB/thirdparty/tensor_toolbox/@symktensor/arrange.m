function X = arrange(X,perm)
%ARRANGE Arranges the rank-1 components of a symktensor.
%
%   ARRANGE(X) normalizes the columns of the factor matrices and then sorts
%   the components by magnitude, greatest to least.
%
%   ARRANGE(X,P) rearranges the components of X according to the
%   permutation P. P should be a permutation of 1 to NCOMPONENTS(X). 
%
%   See also SYMKTENSOR, NCOMPONENTS.
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


%% Just rearrange and return if second argument is a permutation
if exist('perm','var') && isvector(perm)
    X.lambda = X.lambda(perm);
    X.u = X.u(:,perm);
    return;
end

%% Ensure that matrices are normalized
X = normalize(X);

%% Sort
[~, idx] = sort(abs(X.lambda), 1, 'descend');
X.lambda = X.lambda(idx);
X.u = X.u(:,idx);


