function A = subsasgn(A, S, B)
%SUBSASGN Subassignment for symtensor.
%
%   Examples
%   X = symtensor(1:20,3,4);
%   X(:) = (20:-1:1)';  %<- Assign all values
%   X.val = (1:20)'; %<- Another way to assign values
%   X((1:2)') = [-1 -2] %<- Set the first two elements in array
%   X(1,2,3) = -6 %<- Set the value at index (1,2,3) to -6
%   X([4 3 4; 4 4 4]) = -1 * X([4 3 4; 4 4 4]); %<- Reverse the sign
%
%   Note: It is not recommended to assign the same element twice, e.g.,
%   X([5;5]) = [7;8] will result in X(5) = 8. But this behavior is not
%   guaranteed nor tested.
%
%   See also SYMTENSOR, TENSOR/SUBSASGN
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

switch S(1).type
    case '{}'
        error('Cell contents reference from a non-cell array object.')
        
    case '.'
        if isequal(S(1).subs,'val')
            if ~isequal(size(B),size(A.val))
                error('Cannot change the size of ''val'' array');
            end
            A.val = B;
        else
            error('Invalid assignment');
        end
        
    case '()'
        if (A.m > 1) && (numel(S.subs) == A.m)
            if ~isscalar(B)
                error('Can only assign scalars when using subindex');
            end
            if ~all(cellfun(@(x) isscalar(x) && isnumeric(x), S.subs))
                error('Invalid indexing for symktensor');
            end
            newS = S;
            newS.subs = cell(1,1);
            newS.subs{1} = cell2mat(S.subs);
            A = subsasgn(A,newS,B);
        elseif  numel(S.subs) == 1 
            if size(S.subs{1},2) == 1
                A.val = subsasgn(A.val, S, B);
            elseif size(S.subs{1},2) == A.m
                qsubs = sort(S.subs{1},2); %Sort the indices
                asubs = indices(A);
                [~,loca] = ismember(qsubs,asubs,'rows');
                A.val(loca) = B;
            else
                error('Invalid Indexing')
            end
        else
            error('Invalid indexing');
        end        
    otherwise
        %error('Subassignment only allowed for values.');
end

