function [dif, nDif] = aliDif(ali1, ali2)
% Evaluate the difference between two alignment.
%
% Example
%   input   -  ali1.P = [1 1 2 3 4 4 4; ...
%                        1 2 3 4 4 5 6]';
%           -  ali2.P = [1 2 2 2 2 2 3 4 4; ...
%                        1 1 2 3 4 5 5 5 6]';
%   call    -  [dif, nDif] = aliDif(ali1, ali2)
%   output  -  dif = .3, nDif = 8
%
% Input
%   ali1    -  1st alignment
%   ali2    -  2nd alignment
%
% Output
%   dif     -  difference rate
%   nDif    -  difference number between two alignment
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 04-23-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-04-2010

if isstruct(ali1)
    P1 = ali1.P;
    P2 = ali2.P;
else
    P1 = ali1;
    P2 = ali2;
end

n1 = P1(end, 1);
n2 = P1(end, 2);

% Note: this sometimes segfaults. Try running with rowBdSlow to debug.
B1 = rowBd(P1);
B2 = rowBd(P2);

nDif = 0;
for i = 2 : n1
    gap0 = abs(B1(i - 1, 2) - B2(i - 1, 2));
    gap = abs(B1(i, 1) - B2(i, 1));
    
    nDif = nDif + gap0 + .5 * (gap0 ~= gap);
end
nAll = (n1 - 1) * (n2 - 1);

dif = nDif / nAll;
