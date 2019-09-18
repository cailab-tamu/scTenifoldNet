function it = a2it(a)
% Obtain the iteration id from the string.
%
% Input
%   a       -  iteration string, 'tran' | 'corr' | 'wei' | 'mean'
%
% Output
%   it      -  iteration id,        1   |    2   |   3   |   4
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-17-2010

if strcmp(a, 'tran')
    it = 1;

elseif strcmp(a, 'corr')
    it = 2;

elseif strcmp(a, 'wei')
    it = 3;

elseif strcmp(a, 'mean')
    it = 4;

else
    error(['unknown iteration string: ' a]);
end
