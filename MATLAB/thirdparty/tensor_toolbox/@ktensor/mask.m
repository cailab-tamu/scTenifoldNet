function vals = mask(X,W)
%MASK Extract values as specified by a mask tensor.
%
%   V = MASK(X,W) extracts the values in X that correspond to nonzero
%   values in the mask tensor W.
%
%MATLAB Tensor Toolbox.
%Copyright 2017, Sandia Corporation.

% Error check
if any(size(W) > size(X))
    error('Mask cannot be bigger than the input tensor')
end

% Collect info
r = ncomponents(X);
d = ndims(X);
A = X.u; % factor matrices
lambda = X.lambda;

% Extract locations of nonzeros in W
wsubs = find(W);
vsz = [size(wsubs,1) 1];

% Initialize vals array
vals = zeros(vsz);
for j = 1:r
    tmpvals = lambda(j) * ones(vsz);
    for k = 1:d
        akvals = A{k}(wsubs(:,k),j);
        tmpvals = tmpvals .* akvals;
    end
    vals = vals + tmpvals;
end

