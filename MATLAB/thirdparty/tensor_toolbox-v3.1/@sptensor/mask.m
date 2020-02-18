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
    error('Mask cannot be bigger than the data tensor')
end

% Extract locations of nonzeros in W
wsubs = find(W);

% Find which values in the mask match nonzeros in X
[tf,loc] = ismember(wsubs,X.subs,'rows');

% Assemble the final array
nvals = size(wsubs,1);
vals = zeros(nvals,1);
vals(tf) = X.vals(loc(tf));