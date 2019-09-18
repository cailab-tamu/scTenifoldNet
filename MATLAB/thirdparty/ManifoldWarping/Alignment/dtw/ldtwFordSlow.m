function [v, S, DC] = ldtwFordSlow(D,gap_pen)
% Move forward to construct the path matrix.
% Uses Smith-Waterman to do a local alignment
%
% Input
%   D       -  frame (squared) distance matrix, n1 x n2
%
% Output
%   v       -  objective value of dtw
%   S       -  step matrix, n1 x n2
%   DC      -  cumulative distance matrix, n1 x n2
%

assert(gap_pen >= 0, 'Gap penalty must be non-negative');

[n1, n2] = size(D);
S = zeros(n1+1,n2+1);
DC = zeros(n1+1,n2+1);
% init non-start costs
DC(:,1) = (0:n1)*gap_pen;
DC(1,:) = (0:n2)*gap_pen;

for i = 2 : n1+1
    for j = 2 : n2+1
        [DC(i,j), S(i,j)] = min([DC(i,j-1)+gap_pen, DC(i-1,j)+gap_pen, DC(i-1, j-1)]);
        DC(i,j) = DC(i,j) + D(i-1,j-1);
    end
end
%TODO: don't enforce finishing at the end?
v = DC(n1, n2);
