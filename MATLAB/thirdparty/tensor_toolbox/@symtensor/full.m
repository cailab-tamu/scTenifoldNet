function Y = full(X,ver)
%FULL Convert symtensor to a tensor.
%
%   FULL(S) returns a tensor from a symmetric tensor S.
%
%   See also SYMTENSOR, TENSOR.
%
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

%The idx option allows us to pass indexsets as an argument so we don't call
%it each time. Useful if full is called repeatedly on same sized tensors.
%Should probably be deleted because it was broken in the previous release
%anyway (see git repo before 9-1-16)


% Default to new version
if ~exist('ver', 'var');
    ver = 0;
end


switch ver
    case 0 % New version
        
        n = X.n;
        m = X.m;       
        idx = tt_ind2sub(size(X), (1:n^m)');
        classidx = sort(idx, 2);   %Sort indices
        symidx = indices(X);
        [~,refidx] = ismember(classidx, symidx, 'rows');
        newdata = X.val(refidx);
        Y = tensor(reshape(newdata, [size(X) 1]), size(X));
        return;
        
    case 1
        
        I = indices(X);
        sz = X.n * ones(1,X.m);
        Y = tenzeros(sz);
        
        Q = size(I,1);
        
        for q = 1:Q
            i = I(q,:);
            pi = perms(i);
            Y(pi) = X.val(q);
        end
        return;
        
    otherwise;
        error('Incorrect version specification');
end

end






