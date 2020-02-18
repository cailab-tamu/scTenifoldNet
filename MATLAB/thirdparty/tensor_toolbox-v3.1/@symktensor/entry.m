function value = entry(model,I)
%ENTRY Extract a single entry from a symktensor.
%
%   V = ENTRY(M,I) returns the value of M at entry I. M is a symktensor and
%   I is a matrix with rows which are indexes to query.
%   This value is not stored explicitly and so must be computed.
%
%   Examples
%   S = symtenrand(3,4); <-- Declare a random symtensor of size [4,4,4]
%   M = cp_sym(S,2); <-- Decompose S into a symktensor with rank 2
%   I = [1,2,4; 4,4,4];  <-- Matrix of indices to query
%   entry(S,I) <-- Query elements [1,2,4] and [4,4,4]
%
%   See also SYMKTENSOR.
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


%% Extract
lambda = model.lambda;
X = model.u;

%% Case I: I = row
if isrow(I)        
    value = dot(lambda,prod(X(I,:)));
    return;    
end


%% Case II: list of indices
q = size(I,1);
m = size(I,2);
p = size(X,2);
foo = X(I,:);
foo = reshape(foo, [q m p]);
% squeeze(foo(q,:,:)) =  X(I(q,:),:)
bar = squeeze(prod(foo,2));
value = bar * lambda;
