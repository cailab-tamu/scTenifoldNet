function idx = tt_sub2ind64(siz,subs)
%TT_SUB2IND Converts multidimensional subscripts to 64-bit linear indices.
%
%   INDS = TT_SUB2IND(SIZ,SUBS) returns the linear indices
%   equivalent to the subscripts in the array SUBS for a tensor of
%   size SIZ.  
%
%   See also TT_IND2SUB, SUB2IND.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.

if isempty(subs)
    idx = [];
    return;
end

if prod(siz) >= 2^64
    error('Maximum linear index exceeds 2^64');
else
    mult = uint64( [1 cumprod(siz(1:end-1))] );
    idx = uint64( sum( bsxfun( @times, uint64(subs-1), mult), 2 ) + 1 );   
end
