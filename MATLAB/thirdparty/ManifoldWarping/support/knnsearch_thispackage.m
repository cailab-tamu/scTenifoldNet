function idx = knnsearch_thispackage(varargin)
% TODO: try replacing this with the stats toolbox version:
%  - http://www.mathworks.com/help/toolbox/stats/knnsearch.html
%  - also try taking advantage of the case where Q == R 
%
% KNNSEARCH   Linear k-nearest neighbor (KNN) search
% IDX = knnsearch(Q,R,K) searches the reference data set R (d x n array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by each row of Q (d x m array).
% The results are stored in the (m x K) index array, IDX. 
%
% IDX = knnsearch(Q,R) takes the default value K=1.
%
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.

% Check inputs
[Q,R,K,fident] = parseinputs(varargin{:});

% Check outputs
error(nargoutchk(0,2,nargout));

[~,idx] = sort(L2_distance(Q,R));
if fident
    idx = idx(2:K+1,:)';
else
    idx = idx(1:K,:)';
end

function [Q,R,K,fident] = parseinputs(varargin)
% Check input and output
error(nargchk(1,3,nargin));

Q=varargin{1};

if nargin<2
    R=Q;
    fident = true;
else
    fident = false;
    R=varargin{2};
end

if isempty(R)
    fident = true;
    R=Q;
end

if ~fident
    fident = isequal(Q,R);
end

if nargin<3
    K=1;
else
    K=varargin{3};
end
