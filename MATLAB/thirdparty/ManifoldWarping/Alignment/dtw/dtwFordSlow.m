function [v, S, DC] = dtwFordSlow(D)
% Move forward to construct the path matrix.
%
% Example
%   input   -  D = [0 1 1 1 0 1 1 0 1; ...
%                   1 0 1 1 1 1 1 1 0; ...
%                   1 1 0 1 1 1 1 1 1; ...
%                   1 1 1 0 1 0 0 1 1; ...
%                   1 1 1 0 1 0 0 1 1; ...
%                   0 1 1 1 0 1 1 0 1; ...
%                   1 0 1 1 1 1 1 1 0];
%   call    -  [v, P, U] = dtwFord(D)
%   output  -  v = 1;
%              S = [0 1 1 1 1 1 1 1 1; ...
%                   2 3 1 1 1 1 1 1 3; ...
%                   2 2 3 1 1 1 1 1 1; ...
%                   2 2 2 3 1 1 1 1 1; ...
%                   2 2 2 2 1 1 1 1 1; ...
%                   2 2 2 2 3 1 1 3 1; ...
%                   2 2 2 2 2 3 1 2 3]
%              DC = [0 1 2 3 3 4 5 5 6; ...
%                   1 0 1 2 3 4 5 6 5; ...
%                   2 1 0 1 2 3 4 5 6; ...
%                   3 2 1 0 1 1 1 2 3; ...
%                   4 3 2 0 1 1 1 2 3; ...
%                   4 4 3 1 0 1 2 1 2; ...
%                   5 4 4 2 1 1 2 2 1]
%
% Input
%   D       -  frame (squared) distance matrix, n1 x n2
%
% Output
%   v       -  objective value of dtw
%   S       -  step matrix, n1 x n2
%   DC      -  cumulative distance matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-03-2010

[n1, n2] = size(D);
S = zeros(n1,n2);
DC = zeros(n1,n2);

for i = 1 : n1
    for j = 1 : n2
        d = D(i, j);
        if i == 1 && j == 1
            DC(i, j) = d;
            S(i, j) = 0;
        elseif i == 1
            DC(i, j) = DC(i, j - 1) + d;
            S(i, j) = 1;
        elseif j == 1
            DC(i, j) = DC(i - 1, j) + d;
            S(i, j) = 2;
        else
            [DC(i, j), S(i, j)] = min([DC(i, j - 1) + d, DC(i - 1, j) + d, DC(i - 1, j - 1) + d]);
        end
    end
end
v = DC(n1, n2);

