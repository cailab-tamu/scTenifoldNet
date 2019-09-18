function a = it2a(it)
% Obtain the iteration string from the id.
%
% Input
%   it      -  iteration id,        1   |    2   |   3   |   4
%
% Output
%   a       -  iteration string, 'tran' | 'corr' | 'wei' | 'mean'
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-17-2010

as = {'tran', 'corr', 'wei', 'mean'};
a = as{it};
