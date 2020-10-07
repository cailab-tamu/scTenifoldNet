function a = tt_subsubsref(obj,s)
%TT_SUBSUBSREF Helper function for tensor toolbox subsref.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.


if length(s) == 1
    a = obj;
else
    a = subsref(obj, s(2:end));
end

