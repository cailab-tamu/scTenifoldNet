function disp(X, name) 
%DISP Command window display for a symtensor.
%
%   DISP(X) displays a symtensor with no name.
%
%   DISP(X,NAME) displays a tensor with the given name.
%
%   See also SYMTENSOR, SYMTENSOR/DISPLAY.
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

% preallocate
n=X.n;
m=X.m;
sz=length(X.val);
if sz==0   %empty array
   fprintf('%s is an empty symmetric tensor\n', name);
   return
end
output = cell(sz,1);

fprintf('%s is a symmetric tensor with %s modes of dimension %s\n',...
        name, num2str(X.m), num2str(X.n));

I=indices(X);

spc = floor(log10(max(double(I),[],1)))+1;
if numel(spc) == 1
    fmt = ['\t(%' num2str(spc(1)) 'd)%s'];
else
    fmt = ['\t(%' num2str(spc(1)) 'd,'];
    for i = 2:numel(spc)-1
        fmt = [fmt '%' num2str(spc(i)) 'd,'];
    end
    fmt = [fmt '%' num2str(spc(end)) 'd)%s'];
end
%%
% Get values out so that they look nice
savefmt = get(0,'FormatSpacing');
format compact
S = evalc('disp(X.val)');
set(0,'FormatSpacing',savefmt)
S = textscan(S,'%s','delimiter','\n','whitespace','');
S = S{1};
if ~isempty(strfind(S{1},'*'))
    fprintf('%s\n',S{1});
    S = S(2:end);
end
%%
for i = 1:sz
    output{i} = sprintf(fmt, I(i,:) ,S{i});
end
fprintf('%s\n',output{:});

end
    

