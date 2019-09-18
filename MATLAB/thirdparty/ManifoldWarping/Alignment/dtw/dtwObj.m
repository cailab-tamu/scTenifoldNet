function v = dtwObj(Xs, P)
% Measuring the objective of DTW (J_{dtw}).
%
% Input
%   Xs      -  sequences, 1 x 2 (cell), d x ni
%   P       -  alignment path, 2 x m
%
% Output
%   v       -  objective of DTW
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-03-2010

% warping matrix
[W, W1, W2] = P2W(P);

A = Xs{1} * W1' - Xs{2} * W2';
v = trace(A * A');
