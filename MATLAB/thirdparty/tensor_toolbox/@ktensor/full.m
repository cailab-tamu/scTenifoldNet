function t = full(t,vers)
%FULL Convert a ktensor to a (dense) tensor.
%
%   T = FULL(C) converts a ktensor to a (dense) tensor.
%
%   Examples
%   X = ktensor([3; 2], rand(4,2), rand(5,2), rand(3,2));
%   Y = full(A) %<-- equivalent dense tensor
%
%   See also KTENSOR, TENSOR, KTENSOR/DOUBLE.
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

if ~exist('vers','var')
    vers=2;
end

sz = size(t);
if (vers==1)
    data = t.lambda' * khatrirao(t.u,'r')';
else
    % Given d=length(sz) find a partition of modes 1:s and s+1:d that
    % minimizes the memory for the following matrix-multiply.
    minS=minSplit(sz);
    if (minS==length(sz))
        data = khatrirao(t.u,'r') * t.lambda;
    else
        % This unrolls modes 1:minS into rows and minS+1:end into columns
        % of the column-major matrix data which is then converted into a
        % tensor without permutation.
        data = khatrirao(t.u(1:minS),'r') * diag(t.lambda) * khatrirao(t.u(minS+1:end),'r')';
    end
end

t = tensor(data,sz);
end

function [minS]=minSplit(sz)
    % Scan for optimal splitting with minimal memory footprint.
    mLeft=sz(1);
    mRight=prod(sz(2:end));
    minS=1;
    minSum=mLeft+mRight;
    for s=2:length(sz)-1
        mLeft=mLeft*sz(s);
        mRight=mRight/sz(s);
        if (mLeft+mRight<minSum)
            minSum=mLeft+mRight;
            minS=s;
        else
            % Suppose mL >= mR and n >= 1.
            % Then: mL*n + mR/n = mL + mL*(n-1) + mR/n
            %    >= mL + mR*(n-1+1/n) >= mL + mR.
            % 
            % Initially the right term dominates the sum and every factor
            % reduction on the right gives a much smaller increase on the
            % left. Once the left term begins to dominate it will always
            % grow faster than the corresponding reduction on the right.
            break;
        end
    end
end
