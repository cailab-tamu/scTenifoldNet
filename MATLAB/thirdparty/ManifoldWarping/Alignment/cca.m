function [Y1, Y2, V1, V2, me1, me2, v] = cca(X1, X2, par)
% Canonical Correlation Analysis (CCA).
%
% Input
%   X1      -  1st original sample matrix. d1 x n
%   X2      -  2nd original sample matrix. d2 x n
%   par     -  parameter
%     cen   -  flag that the input needs to be centralized, {'y'} | 'n'
%     reg   -  method of regularization, 'fix' | {'cross'} | 'bach'
%              'fix': manually specified
%              'cross': cross-validation
%              'bach': lams is set to [n / 2; n / 2]
%                      See the paper in JMLR 2002, Kernel Independent Component Analysis, for details.
%     lams  -  manually specified regularization value (used if reg == 'fix'), {[0; 0]}
%     k     -  number of k in k-fold cross validation (used if reg == 'cross'), {10}
%     ms    -  number of grid search (used if reg == 'cross'), {[5; 5]}
%     mis   -  minimum value of lambda in search (used if reg == 'cross'), {[1e-4; 1e-4]}
%     mas   -  maximum value of lambda in search (used if reg == 'cross'), {[1e4; 1e4]}
%     sam   -  sampling method to create the grid (used if reg == 'cross'), {'log'} | 'lin'
%     egy   -  percentage of energy to keep for the choice of b, {[]}
%     b     -  num dimensions kept in the reduction
%              if egy == [], then b has be specified.
%
% Output
%   Y1,Y2   -  transformed sample matrices, b x n, b x n
%   V1,V2   -  transformation matrices, d1 x b, d2 x b
%   me1,me2 -  means of sample matrices
%   v       -  objective value of CCA
%

% parameters
if iscell(par), par = cell2option(par); end
if ~isfield(par,'cen'), par.cen = true; end
if ~isfield(par,'reg'), par.reg = 'bach'; end
if ~isfield(par,'lams'), par.lams = [0;0]; end
if strcmp(par.reg, 'cross')
    if ~isfield(par,'k'),   par.k   = 10; end
    if ~isfield(par,'ms'),  par.ms  = [10;10]; end
    if ~isfield(par,'mis'), par.mis = [1e-4;1e-4]; end
    if ~isfield(par,'mas'), par.mas = [1e4;1e4]; end
    if ~isfield(par,'sam'), par.sam = 'log'; end
end
if ~isfield(par,'egy'),  par.egy = []; end
if ~isfield(par,'b'),    par.b = []; end
if ~isfield(par,'debg'), par.debg = false; end

% dimensions
[d1, n] = size(X1);
d2 = size(X2, 1);

% centralize sample matrices
if par.cen
    [X1, me1] = center(X1);
    [X2, me2] = center(X2);
end

% covariance
C11 = X1 * X1';
C22 = X2 * X2';
C12 = X1 * X2';

% regularization
if issparse(C11)
    myrank = @sprank;
else
    myrank = @rank;
end
if myrank(C11) < d1 || myrank(C22) < d2
    switch par.reg
        case 'fix'
            lams = par.lams;
        case 'bach'
            lams = [1; 1] * n / 2;
        case 'cross'
            lams = ccaCross(X1, X2, par);
        otherwise
            error(['unknown regularization algorithm: ' par.reg]);
    end
    if par.debg, fprintf('reg %.5f %.5f\n', lams(1), lams(2)); end;
else
    lams = par.lams;
end

% main algorithm of cca
[V1, V2] = ccaCore(C11, C22, C12, lams, par.egy, par.b);

% projection
Y1 = V1' * X1;
Y2 = V2' * X2;

if nargout > 6
    % objective value
    v = ccaObj(C11, C22, C12, V1, V2, lams);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cX, mu] = center(X)
mu = mean(X, 2);
cX = X - repmat(mu, 1, size(X, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = ccaObj(C11, C22, C12, V1, V2, lams)
% Test usage of CCA.
%
% Input
%   C11   -  covariance, d1 x d1
%   C22   -  covariance, d2 x d2
%   C12   -  covariance, d1 x d2
%   V1    -  1st transformation matrix, d1 x b
%   V2    -  2nd transformation matrix, d2 x b
%   lams  -  regularization weight, 2 x 1
%
% Output
%   Y1    -  1st transformed sample matrix, b x n
%   Y2    -  2nd transformed sample matrix, b x n
%   v     -  objective value of CCA

%[d1, d2] = size(C12);
A = V1' * C12 * V2;
%B1 = V1' * (C11 + lams(1) * eye(d1)) * V1;
%B2 = V2' * (C22 + lams(2) * eye(d2)) * V2;

B1a = V1' * C11 * V1;
B2a = V2' * C22 * V2;

% sneaky way of getting: v = trace((B1a + B2a) \ A)
a = trace(A);
b1 = trace(B1a);
b2 = trace(B2a);
if abs(b1) < 1e-9 || abs(b2) < 1e-9
    v = 0;
else
    v = a / (sqrt(b1) * sqrt(b2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V1, V2] = ccaCore(C11, C22, C12, lams, egy, b)
% Main implementation of CCA.
%
% Input
%   C11   -  covariance, d1 x d1
%   C22   -  covariance, d2 x d2
%   C12   -  covariance, d1 x d2
%   lams  -  regularization weight, 2 x 1
%   egy   -  energy threshold
%   b     -  #dimension to keep
%
% Output
%   V1    -  1st transformation matrix, d1 x b
%   V2    -  2nd transformation matrix, d2 x b

% global
global C11G C22G C12G;
isG = isempty(C11);

% covariance
if isG
    [d1, d2] = size(C12G);
    C1 = [zeros(d1, d1), C12G; ...
        C12G', zeros(d2, d2)];
    C2 = [C11G + lams(1) * eye(d1), zeros(d1, d2); ...
        zeros(d2, d1), C22G + lams(2) * eye(d2)];

else
    %TODO: be smarter about sparsity here?
    [d1, d2] = size(C12);
    C1 = [zeros(d1, d1), C12; ...
        C12', zeros(d2, d2)];
    C2 = [C11 + lams(1) * eye(d1), zeros(d1, d2); ...
        zeros(d2, d1), C22 + lams(2) * eye(d2)];
end

% generalized eigen-decomposition
if issparse(C1)
    [V,D] = eig(full(C1),full(C2)); %eigs(C1, C2, b);
else
    [V,D] = eig(C1, C2);
end
[Lamb, index] = sort(diag(D), 'descend');
V = V(:, index);

% energy
if ~isempty(egy)
    d = min([d1, d2, length(Lamb)]);
    Lamb = Lamb(1 : d);
    b = find(cumsum(Lamb/sum(Lamb)) >= egy, 1);
end

if isempty(b)
    error('must specify either params.b or params.egy for CCA');
end

V1 = V(1 : d1, 1 : b);
V2 = V(d1 + 1 : end, 1 : b);

% when using eigs, sometimes the sign of the eigenvectors flip
% ensure that the first element is positive and everyone is happy
if V1(1)<0, V1 = V1*-1; end
if V2(1)<0, V2 = V2*-1; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lams = ccaCross(X1, X2, par)
% Cross validation for computing the regularization weights.
%
% Input
%   X1    -  1st original sample matrix. d1 x n
%   X2    -  2nd original sample matrix. d2 x n
%  par:
%   k     -  the number of k in k-fold cross validation
%   ms    -  number of grid search if reg == 'cross', 2 x 1
%   mis   -  minimum value of lambda in search, 2 x 1
%   mas   -  maximum value of lambda in search, 2 x 1
%   sam   -  sampling method to create the grid, 'log' | 'lin'
%   egy   -  energy threshold
%   b     -  #dimension to keep
%
% Output
%   lams  -  regularization weights, 2 x 1

% global
global C11G C22G C12G;

% randomly split data points into k fold
n = size(X1,2);
[~, idx] = divN(n, par.k);

% grid of regularization weight
grids = cell(1, 2);
mis = par.mis; mas = par.mas; ms = par.ms;
for i = 1 : 2
    switch par.sam
        case 'log'
        mis(i) = log10(mis(i));
        mas(i) = log10(mas(i));
        grids{i} = logspace(mis(i), mas(i), ms(i));
        case 'lin'
        grids{i} = linspace(mis(i), mas(i), ms(i));
        otherwise
        error(['unknown sampling method: ' par.sam]);
    end
end

% k-fold
Objs = zeros(ms(1), ms(2), par.k);
Bs = zeros(size(Objs));
for c = 1 : par.k
    visTr = idx ~= c;
    visTe = idx == c;

    % sample for training
    X1Tr = X1(:, visTr);
    X2Tr = X2(:, visTr);
    X1Tr = center(X1Tr);
    X2Tr = center(X2Tr);

    % sample for testing
    X1Te = X1(:, visTe);
    X2Te = X2(:, visTe);
    X1Te = center(X1Te);
    X2Te = center(X2Te);

    % covariance
    C11G = X1Tr * X1Tr';
    C22G = X2Tr * X2Tr';
    C12G = X1Tr * X2Tr';
    C11GTe = X1Te * X1Te';
    C22GTe = X2Te * X2Te';
    C12GTe = X1Te * X2Te';

    % grid search for b
    for i1 = 1 : ms(1)
        for i2 = 1 : ms(2)

            % lambda
            lams = [grids{1}(i1); grids{2}(i2)];

            % cca
            V1 = ccaCore([], [], [], lams, par.egy, par.b);

            Bs(i1, i2, c) = size(V1, 2);
        end
    end

    % grid search for obj
    for i1 = 1 : ms(1)
        for i2 = 1 : ms(2)

            % lambda
            lams = [grids{1}(i1); grids{2}(i2)];

            % cca
            [V1, V2] = ccaCore([], [], [], lams, [], min(vec(Bs(:, :, c))));

            % objective value
            Objs(i1, i2, c) = ccaObj(C11GTe, C22GTe, C12GTe, V1, V2, lams);
        end
    end
end
Obj = sum(Objs, 3) / k;

% remove Nan
Obj(isnan(Obj)) = 0;

% debug
if par.debg
    rows = 1; cols = 2;
    axs = iniAx(10, rows, cols, [400 * cols, 400 * rows]);

    [GridX, GridY] = meshgrid(grids{1}, grids{2});

    set(gcf, 'CurrentAxes', axs{1});
    mesh(GridX, GridY, Obj);
    set(axs{1}, 'View', [135 45]);
    title('regularization');

    set(gcf, 'CurrentAxes', axs{2});
    contour(GridX, GridY, Obj);

    if strcmp(sam, 'log')
        set(axs{1}, 'XScale', 'log', 'YScale', 'log');
        set(axs{2}, 'XScale', 'log', 'YScale', 'log');
    else
        set(axs{1}, 'XScale', 'linear', 'YScale', 'linear');
        set(axs{2}, 'XScale', 'linear', 'YScale', 'linear');
    end
end

% optimum value
[~, p] = max(Obj(:));
[p1, p2] = ind2sub(ms, p);
lams = [grids{1}(p1); grids{2}(p2)];
