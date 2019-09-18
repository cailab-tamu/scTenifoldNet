function obj = ctwObj(X1, X2, V1, V2, W1, W2)
% Measuring the objective of CTW.
% \| V1' * X1 * W1 - V2' * X2 * W2 \|_F^2
%
% Input
%   X1      -  1st sequence, d1 x n1
%   X2      -  2nd sequence, d2 x n2
%   V1      -  1st projection, d1 x b
%   V2      -  2nd projection, d2 x b
%   W1      -  1st warping, n1 x n
%   W2      -  2nd warping, n2 x n
%
% Output
%   obj     -  objective of CTW
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-05-2010

if isempty(V1)
    Z1 = X1 * W1;
    Z2 = X2 * W2;
else
    Z1 = V1' * X1 * W1;
    Z2 = V2' * X2 * W2;
end

Z = Z1 - Z2;
obj = sum(Z(:) .^ 2);
