% Check that community components are connected
% Author: Erwan Le Martelot
% Date: 11/11/11
%
% Input
%   - adj : adjacency matrix
%   - coms: communities
%
% Output
%   - nc  : not connected communities

function [nc] = check_coms(adj, coms)
    nc = [];
    for p=1:length(coms)
        ucom = unique(coms{p});
        for c=1:length(ucom)
            idx = coms{p} == c;
            comp = adj(idx,idx);
            if ~is_connected(comp)
                nc = [nc; p c];
                fprintf('Partition %g: Community %g is not connected\n', p, c);
            end
        end
    end
end

