function disp(t,name)
%DISP Command window display of a sumtensor.
%
%   DISP(T) displays a sumtensor with no name.
%
%   DISP(T,NAME) display a sumtensor with the given name.
%
%   See also SUMTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2015, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2015) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt

if ~exist('name','var')
    name = 'ans';
end

if isempty(t.part)
    fprintf(1,'%s is an empty sumtensor\n', name);
    return;
end

fprintf(1,'%s is a sumtensor of size %s with %d parts\n', name, tt_size2str(size(t)), length(t.part));
for i = 1:length(t.part)
    subname = sprintf('%s.part{%d}',name,i);
    disp(t.part{i},subname);
end


