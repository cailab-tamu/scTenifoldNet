% Fast multi-scale detection of communities using Arenas et al.'s parameter
% described in (Arenas, Fernandez, Gomez; Analysis of the structure of
% complex networks at different resolution levels, New Journal of Physics,
% 2008)
% 
% Author: Erwan Le Martelot
% Date: 09/11/11
%
% Input
%   - adj: adjacency matrix
%   - ps : list of r values to consider (from high to low)
%
% Output
%   - coms: community (crip) membership per node
%   - Qs  : criterion values of the corresponding partitions
%
% Comments:
% - This implementation uses Arenas et al. criterion. The code specific to
% it is placed between comment signs lines and can be replaced with code
% optimising any other criterion. It was placed within the algorithm and
% not in external functions for speed optimisation purposes only.

function [coms,Qs] = mscd_afg(adj, ps)

    % Check there is at least one scale parameter
    if (nargin < 2) || isempty(ps)
        error('One scale parameter value at least is required: ms_afg(adj,ps)');
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
        %fprintf('r = %g\n',p);
                
        % Compute current adj
        adj_p = adj + p*eye(length(adj));

        % Total weight m and its double m2
        m2 = sum(sum(adj_p));
        m = m2/2;
        
        % Degree vector
        d = sum(adj_p,2);
        
        % Total weight of a community
        ucom = unique(com);
        for i=1:length(ucom)
            wcom(i) = sum(d(com==ucom(i)));
        end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Initial Q value
        Q = compute_Q(adj_p, com, m2, d);
        
        % While changes can be made
        check_nodes = true;
        check_communities = true;
        while check_nodes
%            disp('AFG big loop');

            % Nodes moving
            moved = true;
            while moved
                %fprintf('Node loop\n');
                moved = false;

                % Create list of nodes to inspect
                l = 1:length(adj_p);

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
                    sum_nb1 = -sum(adj_p(n,nb1));
                    w1 = wcom(com(n)) - d(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for i=1:length(ncom)
                        % Compute dQ for moving n to current community
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        c = ncom(i);
                        nb2 = nb(com(nb) == c);
                        dQ = sum_nb1+sum(adj_p(n,nb2));
                        dQ = (dQ + (d(n)*(w1-wcom(c)))/m2)/m;
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
                        % Update total weight of communities
                        wcom(com(n)) = wcom(com(n)) - d(n);
                        wcom(new_c) = wcom(new_c) + d(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Update community of n
                        com(n) = new_c;
                        % Update Q
                        Q = Q + best_dQ;
                        %%%%%%%%%% DEBUG %%%%%%%%%
% %                         %fprintf('Move node %g: Q=%g\n',n,Q);
%                         eqQ = compute_Q(adj_p, com, m2, d);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    sum_dn1 = sum(d(ncn));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for ncom_idx=1:length(ncom)
                        % Compute dQ for merging cn with current community
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        n2 = com==ncom(ncom_idx);
                        dQ = (sum(sum(adj_p(ncn,n2))) - sum_dn1*sum(d(n2))/m2)/m;
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
                        % Update total weight of communities
                        wcom(new_cn) = wcom(new_cn) + wcom(cn);
                        wcom(cn) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Merge communities
                        com(ncn) = new_cn;
                        % Update Q
                        Q = Q + best_dQ;
                        %%%%%%%%%% DEBUG %%%%%%%%%
% %                         %fprintf('Merge communities %g and %g: Q=%g\n',cn,new_cn,Q);
%                         eqQ = compute_Q(adj_p, com, m2, d);
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
        Qs(p_idx) = Q;

    end % for p

end

% Compute the value of Q for the given partition
function [Q] = compute_Q(adj, com, m2, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = 0;
    for i=1:length(adj)
        Q = Q + adj(i,i);
        for j=i+1:length(adj)
            if com(i) == com(j)
                Q = Q + 2*(adj(i,j) - (d(i)*d(j))/m2);
            end
        end
        Q = Q - (d(i)*d(i))/m2;
    end
    Q = Q / m2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
