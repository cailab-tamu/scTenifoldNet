function B = rowBdSlow(P)
% Bound the correspondence in rows.
%
% Example
%   input   -  P = [1 1 2 3 4 4 4; ...
%                   1 2 3 4 4 5 6]';
%   call    -  B = rowBd(C)
%   output  -  B = [1 3 4 4; ...
%                   2 3 4 6]';
%
% Input
%   P       -  correspondence matrix, n0 x 2
%
% Output
%   B       -  boundary of each row, n1 x 2

n0 = size(P, 1);
n1 = P(end, 1);
B = zeros(n1, 2);

head = 1;
while head <= n0
    i = P(head, 1);
    
    tail = head;
    while tail <= n0 && P(tail, 1) == i
        tail = tail + 1;
    end

    B(i, :) = P([head, tail - 1], 2)';
    
    head = tail;
end
