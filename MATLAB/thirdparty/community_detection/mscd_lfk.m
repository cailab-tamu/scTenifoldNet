% Fast overlapping multi-scale detection of communities using Lancichinetti
% et al.'s criterion described in (Lancichinetti, Fortunato, Kertész;
% Detecting the overlapping and hierarchical community structure in complex
% networks, New Journal of Physics, 2009)
% 
% Author: Erwan Le Martelot
% Date: 11/11/11
%
% Input
%   - adj: adjacency matrix
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
% - This implementation uses Lancichinetti et al's criterion. The code
% specific to it is placed between comment signs lines and can be replaced
% with code optimising any other criterion. It was placed within the
% algorithm and not in external functions for speed optimisation purposes
% only.

function [coms,Qs] = mscd_lfk(adj, ps, merge_th)

    % Check there is at least one scale parameter
    if (nargin < 2) || isempty(ps)
        error('One scale parameter value at least is required: ms_lfk(adj,ps)');
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
    
    % Compute community partition for the current parameter
    pcoms = [];
    for p_idx=1:length(ps)
        
        % Current parameter value
        p = ps(p_idx);

        %fprintf('a = %g\n',p);
        
        % While changes can potentially be made
		grow_coms = true;
        merge_com = true;
		while (grow_coms)
            %fprintf('LFK big loop\n');
            
            % Grow communities
            if isempty(pcoms) % If no comunity exist
                %fprintf('Grow initial communities\n');
                % Initially all nodes are potential seeds
                seeds = 1:length(adj);
                % While there are seed nodes, more communities need to be
                % grown from a random seed
                c = 0;
                while ~isempty(seeds)
                    % Initialise new community
                    idx = randi(length(seeds));
                    c = c + 1;
                    n = seeds(idx);
                    pcoms{c} = n;
                    com_nbs{c} = Nbs{n};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Degrees in out out of the community
                    kin(c) = 0;
                    ktot(c) = sum(adj(n,Nbs{n}));
                    % Grow comunity
                    [pcoms{c}, com_nbs{c}, kin(c), ktot(c), ~] = grow_community(pcoms{c}, com_nbs{c}, kin(c), ktot(c), adj, Nbs, p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Update seeds
                    seeds = setdiff(seeds, pcoms{c});
                    seeds(seeds==n) = [];
                end
            % Otherwise grow existing communities    
            else
                %fprintf('Grow communities\n');
                for c=1:length(pcoms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [pcoms{c}, com_nbs{c}, kin(c), ktot(c), changed] = grow_community(pcoms{c}, com_nbs{c}, kin(c), ktot(c), adj, Nbs, p);
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
                %fprintf('Merge communities\n');
                done = true;
                c1 = 1;
                while c1 <= length(pcoms)
                    c2 = c1 + 1;
                    while c2 <= length(pcoms)
                        com12 = intersect(pcoms{c1}, pcoms{c2});
                        if ~isempty(com12) && ...
                            ((length(com12)/length(pcoms{c1}) >= merge_th) || (length(com12)/length(pcoms{c2}) >= merge_th))
                            pcoms{c1} = union(pcoms{c1}, pcoms{c2});
                            pcoms(c2) = [];
                            com_nbs{c1} = setdiff(union(com_nbs{c1},com_nbs{c2}), pcoms{c1});
                            com_nbs(c2) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            kin(c1)= sum(sum(adj(pcoms{c1},pcoms{c1})));
                            kin(c2) = [];
                            ktot(c1) = sum(sum(adj(pcoms{c1},:)));
                            ktot(c2) = [];
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
            Q = Q + kin(c)/(ktot(c)^p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if ~isempty(pcoms)
            Q = Q/length(pcoms);
        end
        Qs(p_idx) = Q;

    end % for p
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [com, com_nbs, kin, ktot, changed] = grow_community(com, com_nbs, kin, ktot, adj, Nbs, a)
    changed = false;
    done = false;
    while ~done
        done = true;
        % Community fitness
        gfit = kin/(ktot^a);
        % Look for best new fitness with a change
        best_fit = gfit;
        for i=1:length(com_nbs)
            n = com_nbs(i);
            nin = intersect(Nbs{n},com);
            knin = 2*sum(adj(n,nin));
            kntot = sum(adj(n,Nbs{n}));
            fit = (kin+knin)/((ktot+kntot)^a);
            if fit > best_fit
                best_fit = fit;
                best_n = n;
                best_knin = knin;
                best_kntot = kntot;
            end
        end
        % If best fit improves community fitness, add it
        if best_fit > gfit
            %fprintf('BF=%g, F=%g\n', best_fit, gfit);
            com = [com best_n];
            com_nbs(com_nbs==best_n) = [];
            com_nbs = unique([com_nbs setdiff(Nbs{best_n},com)]);
            kin = kin + best_knin;
            ktot = ktot + best_kntot;
            done = false;
            changed = true;
        end
        % If a change has been made to the community, check nodes
        if ~done
            stable = false;
            changed = false;
            while ~stable
                stable = true;
                gfit = kin/(ktot^a);
                for i=1:length(com)
                    n = com(i);
                    knin = 2*sum(adj(n,intersect(Nbs{n},com)));
                    kntot = sum(adj(n,Nbs{n}));
                    fit = (kin-knin)/((ktot-kntot)^a);
                    % If fitness without node is better, remove it
                    if fit > gfit
                        comp = com;
                        comp(comp==n) = [];
                        if is_connected(adj(comp,comp))
                            n = com(i);
                            com(com == n) = [];
                            kin = kin - knin;
                            ktot = ktot - kntot;
                            stable = false;
                            changed = true;
                            % Restart while loop
                            break;
                        end
                    end
                end
            end
            % If some changed occured, recompute the list of neighbours
            if changed
                com_nbs = setdiff(unique([Nbs{com}]), com);
            end
        end
    end % While not done for community
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method inspired from Huang et al.'s work on structural similarity with
% their algorithm applied to Lancichinetti et al's criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [com, com_nbs, kin, ktot, changed] = grow_community2(com, com_nbs, kin, ktot, adj, Nbs, a)
    changed = false;    
    % List of neighbours to explore
    nbs = com_nbs;
    % While more neighbours to explore
    while ~isempty(nbs)
        % Community fitness
        gfit = kin/(ktot^a);
        % Look for best new fitness with a change
        best_fit = gfit;
        for i=1:length(nbs)
            n = nbs(i);
            knin = 2*sum(adj(n,intersect(Nbs{n},com)));
            kntot = sum(adj(n,Nbs{n}));
            fit = (kin+knin)/((ktot+kntot)^a);
            if fit > best_fit
                best_fit = fit;
                best_n = n;
                best_knin = knin;
                best_kntot = kntot;
            end
        end
        % If best fit improves community fitness, add it
        if best_fit > gfit
            %fprintf('BF=%g, F=%g\n', best_fit, gfit);
            com = [com best_n];
            com_nbs(com_nbs==best_n) = [];
            nout = setdiff(Nbs{best_n},com);
            com_nbs = union(com_nbs, nout);
            kin = kin + best_knin;
            ktot = ktot + best_kntot;
            nbs = union(nbs, nout);
            nbs(nbs==best_n) = [];
            changed = true;
        else
            nbs = [];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
