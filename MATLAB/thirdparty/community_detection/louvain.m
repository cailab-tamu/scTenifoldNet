% Louvain method
%
% Input
%   - adj: adjacency matrix
%
% Output
%   - com: communities (listed for each node)
%   - Q  : stability value of the given partition
%
% Author: Erwan Le Martelot
% Date: 01/06/11

function [com, Q, comty] = louvain(adj)
    
    % Run Louvain method
    tic;
    comty = cluster_jl_cpp(adj,1,1,0,0);
    toc;
    
    % Return the partition with highest Q
    [Q, idx] = max(comty.MOD);
    com = comty.COM{idx}';
    
end
