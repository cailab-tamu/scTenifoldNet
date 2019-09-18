% Modularity optimisation based on a greedy agglomerative method
%
% Input
%   - adj: (symmetrical) adjacency matrix
%
% Output
%   - com: communities (listed for each node)
%   - Q  : modularity value of the given partition
%
% Author: Erwan Le Martelot
% Date: 23/11/10

function [com,Q] = fast_newman(adj)

    % Set initial communities with one node per community
    cur_com = [1:length(adj)]';

    % Initialise best community to current value
    com = cur_com;

    % Compute initial community matrix
    e = get_community_matrix(adj,com);
    % Lines and columns sum (speed optimisation)
    ls = sum(e,2);
    cs = sum(e,1);

    % Initialise best known and current Q values
    cur_Q = trace(e) - sum(sum(e^2));
    Q = cur_Q;

    % Loop until no more aggregation is possible
    while length(e) > 1

        % Print progress
        %fprintf('Loop %d/%d...',length(adj)-length(e)+1,length(adj));
        %tic

        % Best Q variation
        loop_best_dQ = -inf;

        % For all the pairs of nodes that could be merged
        can_merge = false;
        for i=1:length(e)
            for j=i+1:length(e)
                % If they share edges
                if e(i,j) > 0
                    % Compute the variation in Q
                    dQ = 2 * (e(i,j) - ls(i)*cs(j));
                    % If best variation, then keep track of the pair
                    if dQ > loop_best_dQ
                        loop_best_dQ = dQ;
                        best_pair = [i,j];
                        can_merge = true;
                    end
                end
            end
        end
        if ~can_merge
            disp('!!! Graph with isolated communities, no more merging possible !!!');
            break;
        end

        % Merge the pair of clusters maximising Q
        best_pair = sort(best_pair);
        for i=1:length(cur_com)
            if cur_com(i) == best_pair(2)
                cur_com(i) = best_pair(1);
            elseif cur_com(i) > best_pair(2)
                cur_com(i) = cur_com(i) - 1;
            end    
        end

        % Update community matrix
        % Slow way (for precision comparison)
        %e = get_community_matrix(adj,com);
        % Faster way
        e(best_pair(1),:) = e(best_pair(1),:) + e(best_pair(2),:);
        e(:,best_pair(1)) = e(:,best_pair(1)) + e(:,best_pair(2));
        e(best_pair(2),:) = [];
        e(:,best_pair(2)) = [];
        % Update lines/colums sum
        ls(best_pair(1)) = ls(best_pair(1)) + ls(best_pair(2));
        cs(best_pair(1)) = cs(best_pair(1)) + cs(best_pair(2));
        ls(best_pair(2)) = [];
        cs(best_pair(2)) = [];

        % Update Q value
        cur_Q = cur_Q + loop_best_dQ;

        % Check consistency
        %eqQ = trace(e) - sum(sum(e*e));
        %if Q ~= eqQ
        %    fprintf('Warning: found Q=%d, should be Q=%d. Diff = %d\n',Q,eqQ,abs(Q-eqQ));
        %end

        % If new Q is better, save current partition
        if cur_Q > Q
            Q = cur_Q;
            com = cur_com;
        end

        %fprintf(' completed in %f(s)\n',toc);

    end

end

