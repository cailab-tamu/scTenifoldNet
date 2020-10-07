function [best_score, A, flag, best_perm] = score(A,B,varargin)
%SCORE Checks if two symktensors match except for permutation.
%   
%   SCORE(A,B) returns the score of the match between A and B where
%   A is trying to be matched against B. It converts both to single-mode
%   ktensors and calls the ktensor SCORE function. 
%
%   See also SYMKTENSOR, KTENSOR/SCORE.
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

%% Make sure A and B are symmetric ktensors
if ~isa(A,'symktensor') || ~isa(B,'symktensor')
    error('Both arguments must be symktensors');
end

A = normalize(A);
B = normalize(B);

AA = ktensor(A.lambda, A.u);
BB = ktensor(B.lambda, B.u);

[best_score, A, flag, best_perm] = score(AA,BB,varargin{:});
