% Kernighan-Lin refinement process based on stability (Markov chain model)
% To use the algorithm with modularity set ts=1
%
% Input
%   - adj: adjacency matrix
%   - com: communities
%   - ts : vector of Markov times to consider
%   - nb_passes: Max number of passes
%
% Output
%   - cur_com: communities after processing
%   - Q  : stability value of the partition after processing
%
% Author: Erwan Le Martelot
% Date: 28/03/11

function [cur_com, Q] = kernighan_lin(adj, com, ts, nb_passes)

    % Set 'true' to display stuff
    verbose = false;

    if nargin < 3
        ts = 1;
    end
    
    if nargin < 4
        nb_passes = inf;
    end

    if ~isempty(find(com==0, 1))
        %adj, com
        error('kernighan_lin: com cannot contain zeros');
    end

    % Compute initial indicator matrix
    H = get_indicator_matrix(adj, com); %diag(ones(1,length(adj)));
    % Nodes degree (line) vector
    d = sum(adj,2)';
    % Transition matrix M
    %M = inv(diag(d)) * adj;
    M = diag(d) \ adj;
    % Markov chain transition probabilities (line) vector
    pi = d / sum(d);
    % Markov chain transition probabilities diagonal matrix
    PI = diag(pi);
    
    % Initial values
    Q = compute_stability(M,H,pi,PI,ts);
    init_Q = Q;
    cur_com = com;
    
    change = true;
    while change
        change = false;
        
        % For each node
        for i=1:length(adj)
            
            % Best Q found in the current loop
            loop_best_Q = Q;
            
            % For each neighbour i that belongs to another cluster
            for j=i+1:length(adj)

                if (adj(i,j) ~= 0) && (cur_com(i) ~= cur_com(j))
                    
                    % Set com(i) as com(j)
                    com_tmp = cur_com;
                    com_tmp(i) = com_tmp(j);
                    k=1;
                    while k <= length(unique(com_tmp))
                        if ~isempty(find(com_tmp==0, 1))
                            %i, j, k, adj, com_tmp
                            error('kernighan_lin: com_tmp cannot contain zeros');
                        end
                        if isempty(find(com_tmp==k, 1))
                            for l=1:length(com_tmp)
                                if com_tmp(l) > k
                                    com_tmp(l) = com_tmp(l) - 1;
                                end
                            end
                        else
                            k = k + 1;
                        end
                    end
                    for s=1:length(ts)
                        H = get_indicator_matrix(adj, com_tmp);
                    end
                    [tmp_Q] = compute_stability(M, H, pi, PI, ts);
                    if tmp_Q > loop_best_Q
                        loop_best_Q = tmp_Q;
                        best_pair = [i j];
                    end
                    
                    % Set com(j) as com(i)
                    com_tmp = cur_com;
                    com_tmp(j) = com_tmp(i);
                    k=1;
                    while k <= length(unique(com_tmp))
                        if isempty(find(com_tmp==k, 1))
                            for l=1:length(com_tmp)
                                if com_tmp(l) > k
                                    com_tmp(l) = com_tmp(l) - 1;
                                end
                            end
                        else
                            k = k + 1;
                        end
                    end
                    for s=1:length(ts)
                        H = get_indicator_matrix(adj, com_tmp);
                    end
                    [tmp_Q] = compute_stability(M,H,pi,PI,ts);
                    if tmp_Q > loop_best_Q
                        loop_best_Q = tmp_Q;
                        best_pair = [j i];
                    end
                    
                end
            end
            
            % If new Q is better, change community of i
            if loop_best_Q > Q
                Q = loop_best_Q;
                cur_com(best_pair(1)) = cur_com(best_pair(2));
                k=1;
                while k <= length(unique(cur_com))
                    if isempty(find(cur_com==k, 1))
                        for l=1:length(cur_com)
                            if cur_com(l) > k
                                cur_com(l) = cur_com(l) - 1;
                            end
                        end
                    else
                        k = k + 1;
                    end
                end
                if nb_passes > 0
                    nb_passes = nb_passes - 1;
                    change = true;
                end
            end

        end

    end
    
    if verbose && (Q > init_Q)
       fprintf('Kernighan-Lin improvement of Q from %g to %g\n', init_Q, Q);
    end

end

% Computation of stability
function [Q,QV] = compute_stability(M,H,pi,PI,ts)

    Q = 1;
    QV = [];
    for s=1:length(ts)
        ds = ts(s) - floor(ts(s));
        if ds == 0
            Rs = H' * ( PI * (M^ts(s)) - pi' * pi ) * H;
        else
            fts = floor(ts(s));
            RsI = H' * ( PI * (M^fts) - pi' * pi ) * H;
            RsF = H' * ( PI * (M^(fts+1)) - pi' * pi ) * H;
            Rs = (1-ds)*RsI + ds*RsF;
        end
        Q_tmp = trace(Rs);
        if (Q_tmp < Q)
            Q = Q_tmp;
        end
        QV(length(QV)+1) = Q_tmp;
    end

end

