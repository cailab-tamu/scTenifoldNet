% Fast multi-scale detection of communities using Ronhovde and Nussinov's
% parameter described in (Ronhovde, Nussinov; Local resolution-limit-free
% Potts model for community detection, Physical Review E, 2010)
% This method only works for binary networks.
% 
% Author: Erwan Le Martelot
% Date: 10/11/11
%
% Input
%   - adj: adjacency matrix (binary nad symmetric)
%   - ps : list of gamma values to consider (from high to low)
%
% Output
%   - coms: community (crip) membership per node
%   - Qs  : criterion values of the corresponding partitions
%
% Comments:
% - This implementation uses Ronhovde and Nussinov's criterion. The code
% specific to it is placed between comment signs lines and can be replaced
% with code optimising any other criterion. It was placed within the
% algorithm and not in external functions for speed optimisation purposes
% only.

function [coms,Qs] = mscd_rn(adj, ps)

    % Check there is at least one scale parameter
    if (nargin < 2) || isempty(ps)
        error('One scale parameter value at least is required: ms_rn(adj,ps)');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Neighbours list for each node
    for i=1:length(adj)
        Nbs{i} = find(adj(i,:));
        Nbs{i}(Nbs{i}==i) = [];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initial community partition: each node in one separate community
    com = 1:length(adj);
 
    % Compute community partition for the current parameter
    for p_idx=1:length(ps)
        
        % Current parameter value
        p = ps(p_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %fprintf('g = %g\n',p);      
        
        % Number of nodes in each community
        ucom = unique(com);
        for i=1:length(ucom)
            scom(i) = sum(com==ucom(i));
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Initial Q value
        Q = compute_Q(adj, com, p);
        
        % While changes can be made
        check_nodes = true;
        check_communities = true;
        while check_nodes
%            disp('RN big loop');

            % Nodes moving
            moved = true;
            while moved
                %fprintf('Node loop\n');
                moved = false;

                % Create list of nodes to inspect
                l = 1:length(adj);

                % While the list of candidates is not finished
                while ~isempty(l)

                    % Pick at random a node n from l and remove it from l
                    idx = randi(length(l));
                    n = l(idx);
                    l(idx) = [];

                    % Find neighbour communities of n
                    ncom = unique(com(Nbs{n}));
                    ncom(ncom == com(n)) = [];
                    % For each neighbour community of n
                    best_dQ = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    nb = Nbs{n};
                    nb1 = nb(com(nb) == com(n));
                    sum_nb1 = sum(adj(n,nb1));
                    sum_nb1 = -sum_nb1 + p*(scom(com(n))-1-sum_nb1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for i=1:length(ncom)
                        % Compute dQ for moving n to current community
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        c = ncom(i);
                        nb2 = nb(com(nb) == c);
                        sum_nb2 = sum(adj(n,nb2));
                        sum_nb2 = sum_nb2 - p*(scom(c)-sum_nb2);
                        dQ = sum_nb1 + sum_nb2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % If best so far, keep track of the best
                        if dQ > best_dQ
                            best_dQ = dQ;
                            new_c = ncom(i);
                        end
                    end

                    % If a move is worth it, do it
                    if best_dQ > 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Update number of nodes in communities
                        scom(com(n)) = scom(com(n)) - 1;
                        scom(new_c) = scom(new_c) + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Update community of n
                        com(n) = new_c;
                        % Update Q
                        Q = Q + best_dQ;
                        %%%%%%%%%% DEBUG %%%%%%%%%
% %                         %fprintf('Move node %g: Q=%g\n',n,Q);
%                         eqQ = compute_Q(adj, com, p);
%                         if abs(Q - eqQ) >= 0.00001
%                             fprintf('Node Warning: found Q=%f, should be Q=%f. Diff = %f\n',Q, eqQ, abs(Q-eqQ));
%                             fprintf('Correcting\n'); Q = eqQ;
%                         end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%
                        % A move occured
                        moved = true;
                        check_communities = true;
                    end

                end

            end % Nodes
            check_nodes = false;
            
            if ~check_communities
                break;
            end
            
            % Community merging
            moved = true;
            while moved
                %fprintf('Community loop\n');
                moved = false;

                % Create community list cl
                cl = unique(com);

                % While the list of candidates is not finished
                while ~isempty(cl)

                    % Pick at random a community cn from cl and remove it from cl
                    idx = randi(length(cl));
                    cn = cl(idx);
                    cl(idx) = []; 

                    % Find neighbour communities of cn
                    ncn = find(com==cn);
                    nbn = unique([Nbs{ncn}]);
                    ncom = unique(com(nbn));
                    ncom(ncom == cn) = [];

                    % For each neighbour community of cn
                    best_dQ = 0;
                    for ncom_idx=1:length(ncom)
                        % Compute dQ for merging cn with current community
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        n2 = com==ncom(ncom_idx);
                        dQ = sum(sum(adj(ncn,n2)));
                        dQ = dQ - p*(scom(cn)*scom(ncom(ncom_idx)) - dQ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % If positive, keep track of the best
                        if dQ > best_dQ
                            best_dQ = dQ;
                            new_cn = ncom(ncom_idx);
                        end
                    end

                    % If a move is worth it, do it
                    if best_dQ > 0
%                        disp('Move com');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Update number of nodes in communities
                        scom(new_cn) = scom(new_cn) + scom(cn);
                        scom(cn) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Merge communities
                        com(ncn) = new_cn;
                        % Update Q
                        Q = Q + best_dQ;
                        %%%%%%%%%% DEBUG %%%%%%%%%
% %                         %fprintf('Merge communities %g and %g: Q=%g\n',cn,new_cn,Q);
%                         eqQ = compute_Q(adj, com, p);
%                         if abs(Q - eqQ) >= 0.00001
%                             fprintf('Com Warning: found Q=%f, should be Q=%f. Diff = %f\n',Q, eqQ, abs(Q-eqQ));
%                             fprintf('Correcting\n'); Q = eqQ;
%                         end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%
                        % A move occured
                        moved = true;
                        check_nodes = true;
                    end

                end

            end % Communities
            check_communities = false;
            
        end % while changes can be made

        % Reindexing communities
        ucom = unique(com);
        for i=1:length(com)
            com(i) = find(ucom==com(i));
        end
        
        % Storing best community found for current parameter value
        coms{p_idx} = com';
        Qs(p_idx) = -Q;

    end % for p

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the value of Q for the given partition
function [Q] = compute_Q(adj, com, g)
    Q = 0;
    ucom = unique(com);
    for i=1:length(ucom)
        ncn = com==ucom(i);
        nn = sum(ncn);
        csum = sum(sum(adj(ncn,ncn)))/2;
        Q = Q + (csum - g*(((nn*(nn-1))/2)-csum));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
