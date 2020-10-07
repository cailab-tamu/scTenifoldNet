function s = tt_intvec2str(v)
%TT_INTVEC2STR Print integer vector to a string with brackets.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.


if isempty(v)
    s = sprintf('[]');
    return;
end

s = ['[ ' sprintf('%d ',v(1:end)) ']'];
