function [P, W] = dtwBackSlow(S)
% Trace back to seek the optimum path.
%
% Example
%   input   -  S = [0 1 1 1 1 1 1 1 1; ...
%                   2 3 1 1 1 1 1 1 3; ...
%                   2 2 3 1 1 1 1 1 1; ...
%                   2 2 2 3 1 1 1 1 1; ...
%                   2 2 2 2 1 1 1 1 1; ...
%                   2 2 2 2 3 1 1 3 1; ...
%                   2 2 2 2 2 3 1 2 3]
%   call    -  [P, W] = dtwBack(S)
%   output  -   P = [1 2 3 4 5 5 5 5 6 7; ...
%                    1 2 3 4 4 5 6 7 8 9]'
%               W = [1 0 0 0 0 0 0 0 0; ...
%                    0 1 0 0 0 0 0 0 0; ...
%                    0 0 1 0 0 0 0 0 0; ...
%                    0 0 0 1 0 0 0 0 0; ...
%                    0 0 0 1 1 1 1 0 0; ...
%                    0 0 0 0 0 0 0 1 0; ...
%                    0 0 0 0 0 0 0 0 1]
%
% Input
%   S       -  step matrix, n1 x n2
%
% Output
%   P       -  warping path matrix, n0 x 2
%   W       -  warping matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-16-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-04-2010

[n1, n2] = size(S);
maM = n1 + n2;
P = zeros(maM, 2);

i = n1; j = n2;
n0 = 1;
P(n0, :) = [i, j];
while true
    s = S(i, j);
    
    if s == 0
        break;
    elseif s == 1
        j = j - 1;
    elseif s == 2
        i = i - 1;
    elseif s == 3
        i = i - 1;
        j = j - 1;     
    else
        error('unknown path direction');
    end

    n0 = n0 + 1;
    P(n0, :) = [i, j];
end
P(n0 + 1 : end, :) = [];
P = P(end : -1 : 1, :);

% warping matrix
Ws = P2W(P);
W = Ws{1}' * Ws{2};
