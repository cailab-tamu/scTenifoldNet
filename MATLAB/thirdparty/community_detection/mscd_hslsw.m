% Fast overlapping multi-scale detection of communities using Huang et al's
% criterion described in (Huang, Sun, Liu, Song, Weninger; Towards Online
% Multiresolution Community Detection in Large-Scale Networks, PLoS ONE,
% 2011)
% 
% Author: Erwan Le Martelot
% Date: 23/11/11
%
% Input
%   - adj: adjacency matrix (symmetrical, can be weighted)
%   - ps : list of alpha values to consider
%   - merge_th: merge threshold (default: 0.5)
%
% Output
%   - coms: list of communities
%   - Qs  : criterion values of the corresponding partitions
%
% Comments:
% - This method uses a local based criterion and hence differs from the
% global criterion approaches in the code, even though the multi-scale
% concept remains the same. Differences are mainly the possibility of
% overlapping communities and the absence of a global criterion replaced
% by local criterion for each cluster.
% - This implementation uses Huang et al's criterion. The code
% specific to it is placed between comment signs lines and can be replaced
% with code optimising any other criterion. It was placed within the
% algorithm and not in external functions for speed optimisation purposes
% only.

function [coms,Qs] = mscd_hslsw(adj, ps, merge_th)

    % Check there is at least one scale parameter
    if (nargin < 2) || isempty(ps)
        error('One scale parameter value at least is required: ms_hslsw(adj,ps)');
    end
    
    % Merge threshold
    if nargin < 3
        merge_th = 0.5;
    end
    
    % Neighbours list for each node
    for i=1:length(adj)
        Nbs{i} = find(adj(i,:));
        Nbs{i}(Nbs{i}==i) = [];
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Similarity between nodes
    for u=1:length(adj)
        for v=u+1:length(adj)
            %Nu = [u Nbs{u}];
            %Nv = [v Nbs{v}];
            if adj(u,v) ~= 0
                INuv = intersect(Nbs{u},Nbs{v});
                %S1(u,v) = length(INuv)/(sqrt(length(Nbs{u})*length(Nbs{v})));
                S(u,v) = ((adj(u,INuv) * adj(v,INuv)') + 2*adj(u,v)) / (sqrt((adj(u,Nbs{u})*adj(u,Nbs{u})')*(adj(v,Nbs{v})*adj(v,Nbs{v})')));
                %if S1(u,v) ~= S(u,v)
                %    fprintf('error\n');
                %end         
                S(v,u) = S(u,v);
            end
        end
    end
%     % Output the similarity matrix
%     for i=1:length(adj)
%         fprintf('%d: ', i);
%         for j=1:length(adj)
%             if adj(i,j) > 0
%                 fprintf('(%d,%f) ', j, S(i,j));
%             end
%         end
%         fprintf('\n');
%     end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute community partition for the current parameter
    pcoms = [];
    for p_idx=1:length(ps)
        
        % Current parameter value
        p = ps(p_idx);
        if debug() fprintf('Current p = %g\n',p); end
        
        % While changes can potentially be made
		grow_coms = true;
        merge_com = true;
		while (grow_coms)
        
            % Grow communities
            if isempty(pcoms) % If no comunity exist
                % Initially all nodes are potential seeds
                seeds = 1:length(adj);
                % While there are seed nodes, more communities need to be
                % grown from a random seed
                c = 0;
                while ~isempty(seeds)
                    if debug() fprintf('Grow initial communities\n'); end
                    % Initialise new community
                    idx = randi(length(seeds));
                    c = c + 1;
                    n = seeds(idx);
                    pcoms{c} = n;
                    com_nbs{c} = Nbs{n};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Structural similarities in out out of the community
                    sin(c) = 0;
                    stot(c) = sum(S(n,Nbs{n}));
                    % Grow comunity
                    if debug() fprintf('Starting community %d from node %d\n', c, n); end
                    [pcoms{c}, com_nbs{c}, sin(c), stot(c), ~] = grow_community(pcoms{c}, com_nbs{c}, sin(c), stot(c), S, Nbs, p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Update seeds
                    seeds = setdiff(seeds, pcoms{c});
                    seeds(seeds==n) = [];
                end
            % Otherwise grow existing communities    
            else
                if debug() fprintf('Grow existing communities (%d)\n', length(pcoms)); end
                for c=1:length(pcoms)
                    if debug() fprintf('Growing community %d\n', c); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [pcoms{c}, com_nbs{c}, sin(c), stot(c), changed] = grow_community(pcoms{c}, com_nbs{c}, sin(c), stot(c), S, Nbs, p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    merge_com = merge_com || changed;
                end
            end
            % Communities need not be grown again unless communites are merged
            grow_coms = false;

            % If no need to merge communities, jump the merging section
			if ~merge_com
                continue;
            end
            
            % Merge 2 communities if they share a significant amount of nodes
            done = false;
            while ~done
                done = true;
                c1 = 1;
                while c1 <= length(pcoms)  
                    c2 = c1 + 1;
                    while c2 <= length(pcoms)
                        com12 = intersect(pcoms{c1}, pcoms{c2});
                        if isempty(com12)
                            c2 = c2 + 1;
                            continue;
                        end
                        if (length(com12)/length(pcoms{c1}) >= merge_th) || (length(com12)/length(pcoms{c2}) >= merge_th)
                            pcoms{c1} = union(pcoms{c1}, pcoms{c2});
                            pcoms(c2) = [];
                            com_nbs{c1} = setdiff(union(com_nbs{c1},com_nbs{c2}), pcoms{c1});
                            com_nbs(c2) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            sin(c1)= sum(sum(S(pcoms{c1},pcoms{c1})));
                            sin(c2) = [];
                            stot(c1) = sum(sum(S(pcoms{c1},:)));
                            stot(c2) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %c1 = length(pcoms)+1; c2 = c1;
                            grow_coms = true;
                            done = false;
                        else
                            c2 = c2 + 1;
                        end
                    end
                    c1 = c1 + 1;
                end
            end
            % Communities need not be checked again unless nodes are moved
			merge_com = false;
            
        end % while grow_com
        
        % Storing communities
        coms{p_idx} = pcoms;
        Q = 0;
        for c=1:length(pcoms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            Q = Q + sin(c)/stot(c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if ~isempty(pcoms)
            Q = Q/length(pcoms);
        end
        Qs(p_idx) = Q;

    end % for p
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [com, com_nbs, sin, stot, changed] = grow_community(com, com_nbs, sin, stot, S, Nbs, a)
    changed = false;
    % List of neighbours to explore
    nbs = com_nbs;
    % While more neighbours to explore
    while ~isempty(nbs)
        % Neighbour with maximum similarity to a node from com
        [~,idx] = max(max(S(nbs,com),[],2));
        n = nbs(idx);
        % Compute its tunable tightness gain
        nin = intersect(Nbs{n},com);
        snin = sum(S(n,nin));
        sntot = sum(S(n,Nbs{n}));
        if sin > 0
            Tn = (stot-sin)/sin - (a*(sntot-snin)-snin)/(2*snin);
        else
            Tn = 1;
        end
        % If tightness gain positive
        if Tn > 0
            if debug() fprintf('Insert %d (Tn=%g > 0)\n', n, Tn); end
            com = [com n];
            com_nbs(com_nbs==n) = [];
            nout = setdiff(Nbs{n}, com);
            com_nbs = union(com_nbs, nout);
            sin = sin + 2*snin;
            stot = stot + sntot;
            nbs = union(nbs, nout);
            changed = true;
        end
        nbs(nbs==n) = [];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Debug switch: on = {true, false}
function [on] = debug()
    on = false;
end
