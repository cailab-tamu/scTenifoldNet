function disp(t, name)
%DISP Command window display for a symktensor.
%
%   DISP(T) displays a symmetric Kruskal tensor with no name.
%
%   DISP(T,NAME) display a symmetric Kruskal tensor with the given name.
%
%   See also DISP, SYMKTENSOR/DISPLAY, SYMKTENSOR
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

fprintf('%s is a symktensor of order %d and dimension %d\n', name, t.m, size(t.u,1));
fprintf('\t%s.lambda = %s\n',name, ['[ ' num2str(t.lambda') ' ]'] );
fprintf('\t%s.U = \n', name);
output = tt_matrix2cellstr(t.u);
fprintf('\t\t%s\n',output{:});
