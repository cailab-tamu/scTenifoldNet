% Create the community matrix from a list giving for each node its community
%
% Input
%   - adj: symmetrical (binary or weighted) adjacency matrix
%   - com: community list
%
% Output
%   - e: community adjacency matrix where e(i,j) is half the sum of the
%        weights of edges connecting communities i and j except for the
%        e(i,i) elements which contain the full sum of the weights of edges
%        connecting community i to itself.
%
% Author: Erwan Le Martelot
% Date: 22/11/10

function [e] = get_community_matrix(adj,com)

    % Number of communities
    nc = length(unique(com));

    % Initialise adjacency matrix list
    e = zeros(nc,nc);

    % Create edges and normalise values by dividing by the sum of all edges
    m = 0;
    for i=1:length(adj)
        for j=i:length(adj)
            if adj(i,j) ~= 0
                ic = com(i);
                jc = com(j);
                if ic == jc
                    e(ic,ic) = e(ic,ic) + adj(i,j);
                else
                    e(ic,jc) = e(ic,jc) + (0.5 * adj(i,j));
                    e(jc,ic) = e(ic,jc);
                end
                m = m + adj(i,j);
            end
        end
    end
    e = e / m;

end

