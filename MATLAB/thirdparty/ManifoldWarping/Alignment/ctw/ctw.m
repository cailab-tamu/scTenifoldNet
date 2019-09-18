function [ali, Ys, Vs, objs, its] = ctw(Xs, ali0, parCtw, parCca)
% Canonical Time Warping (CTW) and Kernel CTW (KCTW)
%
% Input
%   Xs       -  sequences (or similarity matrices), 1 x 2 (cell)
%               If Xs contains sequence, dimension of each Xi is di x ni.
%               If Xs contains similarity, dimension of each Xi is ni x ni.
%   ali0     -  initial alignment
%   parCtw   -  parameter for CTW
%     knl    -  kernel flag, 'y' | {'n'}
%     th     -  stop threshold, {.01}
%     debg   -  debug flag, 'y' | {'n'}
%     nItMa  -  maximum iteration number, {100}
%     PT     -  ground-truth alignment, {[]}
%   parCca   -  parameter for CCA
%               See function cca for more details about the setting of parCca.
%
% Output
%   ali      -  new alignment
%   Ys       -  new sequences, 1 x 2 (cell)
%   Vs       -  transformation, 1 x 2 (cell)
%               If Xs contains sequence, dimension of each Vi is di x b
%               If Xs contains similarity, dimension of each Vi is ni x b
%   objs     -  objective value, 1 x nIt
%   its      -  iteration step ids, 1 x nIt
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 05-20-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 09-06-2010
 
% parameter
if iscell(parCtw), parCtw = cell2option(parCtw); end
isKnl = psY(parCtw, 'knl', 'n');
th = ps(parCtw, 'th', .01);
isDebg = psY(parCtw, 'debg', 'n');
isVdo = psY(parCtw, 'vdo', 'n');
nItMa = ps(parCtw, 'nItMa', 100);
PT = ps(parCtw, 'PT', []);

% debug
if isDebg
    rows = 2; cols = 3;
    fig = 1;
    axs = iniAx(fig, rows, cols, [300 * cols, 300 * rows]);

    if isKnl && isVdo
        vout = vWriter('kctw.avi', 1);
    elseif isVdo
        vout = vWriter('ctw.avi', 1);
    end
end

% initial sequence
X10 = Xs{1};
X20 = Xs{2};

% initial alignment
P0 = ali0.P;
Ws = P2W(P0);
W1 = Ws{1};
W2 = Ws{2};

% coordinate-descent search
objs = zeros(1, nItMa);
its  = zeros(1, nItMa);
nIt = 0;
V1 = []; V2 = [];
while nIt < nItMa
    if isDebg, fprintf('%d\n', nIt); end

    % spatial transformation
    % before transformation, the sample has to be centralized
    if isKnl
        X1 = cenK(X10, W1);
        X2 = cenK(X20, W2);
    else
        X1 = cenX(X10, W1);
        X2 = cenX(X20, W2);
    end
    
    Z1 = X1 * W1;
    Z2 = X2 * W2;
    [~, ~, V1, V2] = cca(Z1, Z2, parCca);
    b = size(V1, 2);

    nIt = nIt + 1;
    objs(nIt) = ctwObj(X1, X2, V1, V2, W1, W2);
    its(nIt) = a2it('tran');
    if isDebg
        debg(axs, objs(1 : nIt), its(1 : nIt), X1, X2, V1, V2, P0, PT);
        if isVdo
            vout = vAddframe(vout);
        end
    end

    % temporal warping
    Y1 = V1' * X1;
    Y2 = V2' * X2;
    D = L2_distance(Y1, Y2);
    [v, S] = dtwFord(D);
    P = dtwBack(S);
    Ws = P2W(P);
    W1 = Ws{1};
    W2 = Ws{2};
    nIt = nIt + 1;  
    objs(nIt) = ctwObj(X1, X2, V1, V2, W1, W2);
    its(nIt) = a2it('corr');
    if isDebg
        debg(axs, objs(1 : nIt), its(1 : nIt), X1, X2, V1, V2, P, PT);
        if isVdo
            vout = vAddframe(vout);
        end
    end

    % stop condition
    dif = aliDif(P, P0);
    if dif <= th
        break;
    end
    P0 = P;
end
objs(nIt + 1 : end) = [];
its(nIt + 1 : end) = [];
if isDebg
    if isKnl, fprintf('k'); end;  % kctw
    fprintf('ctw converged in %d steps\n', nIt);
    if isVdo
        vClose(vout);
    end
end

ali = newAli('P', P, 'obj', objs(end));
if ~isempty(PT)
    ali.dif = aliDif(P, PT);
end
Ys = {Y1, Y2};
Vs = {V1, V2};

%%%%%%%%%%%%%%%%%%%%%%%%
function X = cenX(X0, W)
% Centralize samples.
%
% Input
%   X0  -  original sample matrix, d x nx
%   W   -  warping matrix, nx x n
%
% Output
%   X   -  new sample matrix, d x nx

[nx, n] = size(W);

a = sum(W(:));
b = X0 * W * ones(n, 1);
me = b / a;
X = X0 - repmat(me, 1, nx);


%%%%%%%%%%%%%%%%%%%%%%%%
function K = cenK(K0, W)
% Centralize kernel matrix.
%
% Input
%   K0  -  original similarity matrix, nk x nk
%   W   -  weighted matrix, nk x n
%
% Output
%   K   -  new similarity matrix, nk x nk

[nk, n] = size(W);
a = sum(W(:));

P = eye(nk) - ones(nk, n) * W' / a;
K = P * K0 * P';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function debg(axs, objs, its, X1, X2, V1, V2, P, PT)
% Show current status in the optimization of CTW.

% obj
shIter(objs, its, 'ax', axs{1, 1}, 'mkSiz', 7, 'mkEg', 'n', 'lnWid', 1);
axis square;
if length(objs) > 2
    ylim([0 objs(2) * 1.5]);
end
title('objective');

shSeq({X1, X2}, 'ax', axs{1, 2}, 'mkSiz', 3, 'mkEg', 'n', 'lnWid', 1);
xlabel('x'); ylabel('y');
title('original sequence');

if isempty(V1)
    Z1 = X1;
    Z2 = X2;
else
    Z1 = V1' * X1;
    Z2 = V2' * X2;
end
shSeq({Z1, Z2}, 'ax', axs{1, 3}, 'mkSiz', 3, 'mkEg', 'n', 'lnWid', 1);
xlabel('x'); ylabel('y');
title('new sequence');

shAliChan({X1, X2}, [], 'ax', axs{2, 2}, 'mkSiz', 3, 'mkEg', 'n', 'lnWid', 1, 'nms', {'x', 'y'});
title('original sequence');

ali.P = P;
shAliChan({Z1, Z2}, ali, 'ax', axs{2, 3}, 'mkSiz', 3, 'mkEg', 'n', 'lnWid', 1, 'nms', {'x', 'y'});
title('new sequence');

D = L2_distance(Z1, Z2);
shM(D, 'dis', 'imagesc', 'P', P, 'ax', axs{2, 1}, 'clMap', 'hsv', 'lnMk', '--', 'lnWid', 2, 'lnCl', 'b');
if ~isempty(PT)
    plot(PT(:, 2), PT(:, 1), '-', 'Color', 'k', 'LineWidth', 2);
end
title('Distance');

drawnow;
