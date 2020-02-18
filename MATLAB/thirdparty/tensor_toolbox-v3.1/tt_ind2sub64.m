function subs = tt_ind2sub64(siz,idx)
%TT_IND2SUB Multiple subscripts from linear indices.
%
%   SUBS = TT_IND2SUB(SIZ,INDS) returns that subscripts equivalent
%   to the linear indices in INDS for a tensor of size SIZ.
%
%   See also TT_SUB2IND, IND2SUB, STRATIFIED_SAMPLE.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.

% Created by Tamara G. Kolda, Fall 2018. 


if isempty(idx)
    subs = [];
    return;
end

k = uint64([1 cumprod(siz(1:end-1))]);
n = length(siz);
idx=idx-1;
for i = n : -1 : 1 
    div=idivide(idx,k(i),'floor');
    subs(:,i) = div+1;
    idx=idx-k(i)*div;
end
subs = double(subs);