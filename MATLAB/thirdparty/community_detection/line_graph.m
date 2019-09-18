% Compute the line graph of the given graph
%
% Input
%   - adj: adjacency matrix of the graph
%
% Output
%   - lg_adj  : adjacency matrix of the line graph
%   - lg_nodes: nodes of the initial graph forming each line graph node
%
% Author: Erwan Le Martelot
% Date: 18/02/11

function [lg_adj,lg_nodes] = line_graph(adj)
    lg_nodes = [];
    for i=1:length(adj)
        for j=i+1:length(adj)
            if adj(i,j) > 0
                lg_nodes = [lg_nodes; i j];
            end
        end
    end
    lg_adj = zeros(size(lg_nodes,1),size(lg_nodes,1));
    for i=1:size(lg_nodes,1)
        for j=i+1:size(lg_nodes,1)
            if lg_nodes(i,1) == lg_nodes(j,1) | lg_nodes(i,1) == lg_nodes(j,2) | ...
               lg_nodes(i,2) == lg_nodes(j,1) | lg_nodes(i,2) == lg_nodes(j,2)
                lg_adj(i,j) = 1;
                lg_adj(j,i) = 1;
            end
        end
    end
end
