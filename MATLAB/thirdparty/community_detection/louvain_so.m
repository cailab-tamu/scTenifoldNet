% Stability optimisation based on the Louvain method.
%
% Input
%   - adj: (symmetrical) adjacency matrix
%   - t  : Markov time to consider
%
% Output
%   - com: communities (listed for each node)
%   - Q  : stability value of the given partition
%
% Author: Erwan Le Martelot
% Date: 06/12/10

function [com, Q] = louvain_so(adj, t, model)

    if nargin < 3
        model = 'discrete';
    end

    % Check connections
    single = find(sum(adj)==0);
    
    if ~isempty(single)
        fprintf('Found single nodes: '); disp(single);
        not_single = setdiff(1:length(adj), single);
        adj = adj(not_single,not_single);
    end

    % Compute initial indicator matrix
    H = diag(ones(1,length(adj)));

    % Nodes degree (line) vector
    d = sum(adj,2)';

    % Transition matrix M
    %M = inv(diag(d)) * adj;
    D = diag(d);
    M = D \ adj;
    
    % Markov chain transition probabilities (line) vector
    pi = d / sum(d);
    
    % Compute At
    if strcmp(model,'discrete')
        dt = t - floor(t);
        if dt == 0
            At = D * (M^t);
        else
            ft = floor(t);
            At = D * ((1-dt)*(M^ft) + dt*(M^(ft+1)));
        end
    elseif strcmp(model,'continuous')
        At = D * expm((M-eye(length(M)))*t);
    else
        error('Unknown model');
    end
    
%     % Remove edges that did not exist in A
%     for i=1:length(adj)
%         for j=i+1:length(adj)
%             if adj(i,j) == 0
%                 At(i,j) = 0;
%                 At(j,i) = 0;
%             end
%         end
%     end
    
    % Apply community detection
    [com, Q] = louvain(At);
    
    % Add single nodes back with community 0
    if ~isempty(single)
        tmp = [not_single' com; single' zeros(length(single),1)];
        tmp = sortrows(tmp,1);
        com = tmp(:,2);
    end
    
end