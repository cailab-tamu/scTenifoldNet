function [ali,Ys] = dtw(Xs,D)
% Dynamic Time Warping (DTW).
%
% Input
%   Xs      -  sequences, 1 x 2 (cell), dim x ni
%   D       -  distance matrix (optional, will compute L2 if omitted)
%
% Output
%   ali     -  alignment
%     P     -  warping path, n0 x 2
%   Ys      -  warped inputs
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-09-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-04-2010

if (nargin == 1)
    D = L2_distance(Xs{1}, Xs{2});
end

[v, S] = dtwFord(D);
P = dtwBack(S);

ali = newAli('P', P, 'obj', v);

if nargout > 1
    p = P(:,1);
    q = P(:,2);
    X1 = Xs{1};

    warp_idxs = zeros(1, size(X1,2));
    for i = 1:length(warp_idxs)
        warp_idxs(i) = q(find(p>=i,1));
    end
    warp_idxs = uint32(warp_idxs);

    Ys = {X1, X1(:, warp_idxs)};
end
