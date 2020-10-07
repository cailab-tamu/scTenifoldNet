function U = tocell(X,N)
%TOCELL Convert X to a cell array.
%
%   TOCELL(X) converts X to a cell array, evenly distributing the
%   weight in lambda.
%
%   TOCELL(X,N) absorbs the weights into the Nth factor matrix.
%
%   Examples
%   K = ktensor([2; 4], ones(3,2), ones(4,2), ones(2,2));
%   tocell(K) %<--Output in a cell array with 3 matrices
%   tocell(K,3) %<--Same as above, but weight absorbed into 3rd matrix
%
%   See also KTENSOR, NORMALIZE.
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

if exist('N','var')
    X = normalize(X,N);
    U = X.u;
    return;
end

if isequal(X.lambda,ones(size(X.lambda)))
    U = X.u;
    return;
end

lsgn = sign(X.lambda);
lsplit = nthroot(abs(X.lambda),ndims(X));
U = X.u;
U{1} = U{1} * diag(lsgn);
D = diag(lsplit);
for n = 1:ndims(X)
    U{n} = U{n} * D;
end
        

