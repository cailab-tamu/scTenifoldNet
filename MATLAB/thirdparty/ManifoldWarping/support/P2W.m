function Ws = P2W(P)
% Convert warping path vector to warping path matrix.
%
% Example
%   input   -  P = [1 1; ...
%                   2 2; ...
%                   2 3; ...
%                   3 3; ...
%                   4 3; ...
%                   5 4]
%   call    -  Ws = P2W(P)
%   output  -  Ws{1} = [1 0 0 0 0 0; ...
%                       0 1 1 0 0 0; ...
%                       0 0 0 1 0 0; ...
%                       0 0 0 0 1 0; ...
%                       0 0 0 0 0 1]
%              Ws{2} = [1 0 0 0 0 0; ...
%                       0 1 0 0 0 0; ...
%                       0 0 1 1 1 0; ...
%                       0 0 0 0 0 1]
%
% Input
%   P       -  warping path vector, n x m
%
% Output
%   Ws      -  warping path matrix, 1 x m (cell), ni x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 04-02-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-05-2010

[n, m] = size(P);
ns = P(n, :);

Ws = cell(1, m);
for i = 1 : m
    ni = ns(i);
    Ws{i} = sparse(ni, n);
    Ws{i}(sub2ind([ni, n], P(:, i)', 1 : n)) = 1;
end

if m == 1
    Ws = Ws{1};
end
