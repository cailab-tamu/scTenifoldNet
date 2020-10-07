%% Alternating least squares for Tucker model 
% The function |tucker_als| computes the best rank(R1,R2,..,Rn)
% approximation of tensor X, according to the specified dimensions in
% vector R.  The input X can be a tensor, sptensor, ktensor, or
% ttensor.  The result returned in T is a ttensor.
%
% The method is originally from Tucker (1966) and later revisited in 
% De Lathauwer et al. (2000).
%
% * Tucker, L. R. 
%   Some mathematical notes on three-mode factor analysis.
%   Psychometrika, 1966, 31, 279-311.
% * De Lathauwer, L.; De Moor, B. & Vandewalle, J. 
%   On the best rank-1 and rank-(R_1, R_2, R_N) approximation of
%   higher-order tensors.
%   SIAM Journal on Matrix Analysis and Applications, 2000, 21, 1324-1342.
%
% Note: Oftentimes it's better to use |hosvd| instead.

%% Create a data tensor of size [5 4 3]
rng('default'); rng(0); %<-- Set seed for reproducibility
X = sptenrand([5 4 3], 10)
%% Create a [2 2 2] approximation
T = tucker_als(X,2)        %<-- best rank(2,2,2) approximation 
%% Create a [2 2 1] approximation
T = tucker_als(X,[2 2 1])  %<-- best rank(2,2,1) approximation 
%% Use a different ordering of the dimensions
T = tucker_als(X,2,struct('dimorder',[3 2 1]))
%% Use the n-vecs initialization method
% This initialization is more expensive but generally works very well.
T = tucker_als(X,2,struct('dimorder',[3 2 1],'init','eigs'))
%% Specify the initial guess manually
U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of T
T = tucker_als(X,2,struct('dimorder',[3 2 1],'init',{U0}))
