% Create the indicator matrix from a community list
% giving for each node its community.
%
% Input
%   - adj: symmetrical (binary or weighted) adjacency matrix
%   - com: community list
%
% Output
%   - H: community indicator matrix where H(i,j) is 1 if node i belongs to
%        community j
%
% Author: Erwan Le Martelot
% Date: 06/12/10

function [H] = get_indicator_matrix(adj,com)

    % Number of nodes
    n = length(adj);

    % Number of communities
    nc = length(unique(com));

    % Initialise indicator matrix
    H = zeros(n,nc);

    % Fill in matrix
    for i=1:n
        H(i,com(i)) = 1;
    end

end

