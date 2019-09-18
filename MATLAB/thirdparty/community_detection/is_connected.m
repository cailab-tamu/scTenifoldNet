% Test if a graph is connected
% Author: Erwan Le Martelot
% Date: 11/11/11
%
% Input
%   - adj: adjacency matrix
%
% Output
%   - result: true if connected, false otherwise

function [result] = is_connected(adj)
    % Initially no node visited 
    visited = false(length(adj),1);
    % Visit the first node
    queue = 1;
    visited(1) = true;
    % While there are nodes to visit
    while ~isempty(queue)
        % Take next in queue
        n = queue(1);
        queue(1) = [];
        % For each neighbour of n enqueue them if not visited
        nbs = find(adj(n,:));
        for i=1:length(nbs)
            if ~visited(nbs(i))
                visited(nbs(i)) = true;
                queue = [queue nbs(i)];
            end
        end
    end
    result = sum(visited) == length(adj);
end
