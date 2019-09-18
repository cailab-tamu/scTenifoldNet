% Compute the modularity of a given partition (using stability)
%
% Input
%   - adj: adjacency matrix
%   - com: community list (per node)
%
% Output
%   - Q: Modularity value
%
% Author: Erwan Le Martelot
% Date: 07/04/11

function [Q] = modularity(adj, com)

    Q = stability(adj, com, 1);

end

