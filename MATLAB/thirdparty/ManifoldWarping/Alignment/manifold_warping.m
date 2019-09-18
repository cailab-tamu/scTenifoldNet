function [P,Y1,Y2,V1,V2] = manifold_warping(X1, X2, mode, target_dim, k, mu, thresh, max_its)
% Manifold alignment + DTW
%
% Input
%   X1,X2     - sequences,  each Xi is di x ni.
%   mode      - one of {'linear','nonlinear','embed'}
%
% Output
%   P        -  new alignment path
%   Y1,Y2    -  new sequences
%   V1,V2    -  new mappings (only if mode = 'linear')

% parameters
dim = min(size(X1,1),size(X2,1));
epsilon = 1e-6;
if ~exist('target_dim','var'), target_dim = dim-1; end
if ~exist('k','var'),          k = 10; end
if ~exist('mu','var'),         mu = 1; end
if ~exist('thresh','var'),     thresh = 0.01; end
if ~exist('max_its','var'),    max_its = 100; end

% initial alignment
W1 = createKnnGraph(knnsearch_same(X1, k));
W2 = createKnnGraph(knnsearch_same(X2, k));

if strcmp(mode,'embed')
    X1 = laplacian_eigen(X1', target_dim, k)';
    X2 = laplacian_eigen(X2', target_dim, k)';
    %plot_correspondences(X1,X2);
    if (size(X1,2) ~= size(W1,2) || size(X2,2) ~= size(W2,2))
        size(X1)
        size(X2)
        [size(W1,2) size(W2,2)]
        error('embedding didn''t preserve all instances');
    end
end

is_linear = ~strcmp(mode,'nonlinear');

W12 = sparse(size(X1,2),size(X2,2));
W12(1,1) = 1;                         %Align the 2 first points
W12(size(X1,2),size(X2,2)) = 1;       %Align the 2 last points
%[W12,P0] = my_dtw(X1,X2);

% coordinate-descent search
nIt = 0;
P0 = [(1:size(X1,2))' ones(size(X1,2),1)];  % silly initial alignment
while nIt < max_its
    nIt = nIt + 1;

    % manifold alignment
    if is_linear
        [V1,V2] = manifold_linear(X1,X2,W1,W2,W12,mu,dim+4,epsilon);
        V1 = V1(:,1:target_dim);
        V2 = V2(:,1:target_dim);
        Y1 = V1' * X1;
        Y2 = V2' * X2;
    else
        [Y1,Y2] = manifold_nonlinear(W1,W2,W12,mu,dim+4,epsilon);
        Y1 = Y1(:,1:target_dim)';
        Y2 = Y2(:,1:target_dim)';
    end

    % temporal warping
    [W12,P] = my_dtw(Y1,Y2);

    % stop condition
    dif = aliDif(P, P0);
    fprintf('Iteration %d, warping path change = %g\n',nIt,dif);
    if dif <= thresh
        break;
    end
    P0 = P;
end

if ~is_linear
    V1 = []; V2 = [];
end

end


function [W12,P] = my_dtw(Y1,Y2)
    D = L2_distance(Y1, Y2);
    [~, S] = dtwFord(D);
    P = dtwBack(S);
    ns = P(size(P,1), :);
    W12 = sparse(ns(1),ns(2));
    W12(sub2ind(ns, P(:,1)', P(:,2)')) = 1;
end
