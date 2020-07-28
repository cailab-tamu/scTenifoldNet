function Rpath=FindRpath
% This function finds the path for the installed R (Project for Statistical Computing) in Windows environment.
% Only works in Windows environments.
% e.g.:
% Rpath=FindPathR
% >> 'C:\Program Files\R\R-3.1.1\bin'
if isunix
    [a,b]=system('which R');
    if a==0
    Rpath=deblank(b);
    else
        Rpath='';
    end
return;
end

sep = filesep; env=myGetEnv; 
a=FindWhich(env(:,1),'ProgramFiles');
b=FindWhich(env(:,1),'ProgramW6432'); 
c=FindWhich(env(:,1),'ProgramFiles(x86)');
n=[a;b;c]; isfound=0;
 for i=1:length(n)
     programPath=env{n(i),2};
     D=dir([programPath filesep 'R']);
     if ~isempty(D)
         A={D.name}; B=find(cell2mat(cellfun(@(s) contains(s,'R-'),A,'uniformoutput',0)),1);
         if ~isempty(B),  isfound=1;Rpath=[programPath sep 'R' sep A{B} sep 'bin'];break; end
     end%
 end
if isfound==0,  Rpath=[];end
end % end of FindRpath
function env=myGetEnv
% Get the environmental values from the operation system.
% Returns key-value pairs stored in 'env' cell array.
% e.g.:
%  >> env=myGetEnv;
% >>  env{8,1} =  'ProgramFiles'
%  >> env{8,2} =  'C:\Program Files'
% Weirong Chen  March-8-2015
envmap = java.lang.System.getenv();
envkey=cell(envmap.keySet.toArray);
envvalue=cell(envmap.values.toArray);
env=[envkey envvalue];
end % end of function

function n=FindWhich(StringCellArray, TargetString, nOutput)
% This funcion finds the 'TargetString' in a string cell array and returns
% the serial number of the element found in 'StringCellArray'.
% 'nOutput' = number of output elements.
% Weirong Chen    SEP-14-2013

sn=1:length(StringCellArray); A=cellismember(StringCellArray,TargetString); n=sn(A);
if nargin<3, return;end
if numel(n)>0, n=n(1:nOutput);end
end % FindWhich

function Lia=cellismember(A,B) 
% The built-in "ismember" function in MATLAB fails to perform when the input variables are cells containing different types of variables.
% This function 'cellismember' is a function that performs 'ismember' on
%  cells with various data types.
%  The input A and B must be cell arrays.
%  Example: 
%  Input: A = {'ab','cd', NaN, [], 5, 1}; B = {[], 'cd', NaN, 1};
%  output: Lia = [0 1 1 1 0 1]; 
%
% Acknowledgement: 
% This function greatly benefits from Jan Simon's comments. The previous version was errorful. 
% See 'ismember' for more information
% Weirong Chen   Apr-16-2
% Update: Jun-1-2015015

if ~iscell(B) && ischar(B), B={B}; end% Convert single value into 1x1 cell array of string
str_Index_A = cellfun('isclass', A, 'char');
str_Index_B = cellfun('isclass', B, 'char');
NaN_Index_A=logical(cell2mat(cellfun(@sum, cellfun(@isnan, A,'UniformOutput',false),'UniformOutput',false)));
NaN_Index_B=logical(cell2mat(cellfun(@sum, cellfun(@isnan, B,'UniformOutput',false),'UniformOutput',false)));
empty_Index_A=cellfun(@isempty, A);
empty_Index_B=cellfun(@isempty, B);
num_Index_A=cellfun(@isnumeric, A) & ~NaN_Index_A & ~empty_Index_A; % 'isnumeric' includes NaN and EMPTY.
num_Index_B=cellfun(@isnumeric, B) & ~NaN_Index_B & ~empty_Index_B; % 'isnumeric' includes NaN and EMPTY.
if sum(str_Index_A)>0 && sum(str_Index_B)>0 
    out_Index_str = str_Index_A;
    out_Index_str(str_Index_A)=ismember(A(str_Index_A),B(str_Index_B));
else 
    out_Index_str = false(size(A,1),size(A,2));
end
if sum(num_Index_A)>0 && sum(num_Index_B)>0 
    out_Index_num = num_Index_A;
    out_Index_num(num_Index_A)=ismember(cell2mat(A(num_Index_A)),cell2mat(B(num_Index_B)));
else 
    out_Index_num = false(size(A,1),size(A,2));
end
if sum(NaN_Index_A)>0 && sum(NaN_Index_B)>0 
    out_Index_NaN = NaN_Index_A;
else 
    out_Index_NaN = false(size(A,1),size(A,2));
end
if sum(empty_Index_A)>0 && sum(empty_Index_B)>0 
    out_Index_empty = empty_Index_A;
else 
    out_Index_empty = false(size(A,1),size(A,2));
end
Lia = out_Index_str  | out_Index_num | out_Index_NaN | out_Index_empty;
end %end of cellismember

