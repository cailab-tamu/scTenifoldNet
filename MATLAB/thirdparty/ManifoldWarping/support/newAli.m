function ali = newAli(varargin)
% Create alignment structure.
%
% Input
%   varargin
%     alg   -  algorithm name, {[]}
%     P     -  alignment path, {[]}
%     acc   -  accuracy, {[]}
%     tim   -  time cost, {[]}
%     obj   -  objective value, {[]}
%     dif   -  difference rate (from ground truth alignment), {[]}
%
% Output
%   ali     -  alignment struct
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-09-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-06-2010

alg = ps(varargin, 'alg', []);
P = ps(varargin, 'P', []);
acc = ps(varargin, 'acc', []);
tim = ps(varargin, 'tim', []);
obj = ps(varargin, 'obj', []);
dif = ps(varargin, 'dif', []);

ali = struct('alg', alg, 'P', P, 'acc', acc, 'tim', tim, 'obj', obj, 'dif', dif);
