% Greedy stabilisation optimisation given a model and heuristics
% with possibility of input pre-processing and output post-processing.
% This function is an interface enabling GSO, randomised GSO and multistep
% GSO (i.e. calling gso, rgso and msgso on both discrete and continuous
% models).
%
% Input
%   - adj  : adjacency matrix
%   - ts   : vector of Markov times to consider
%   - model: 'discrete' or 'continuous' Markov model, their randomised 
%            version 'rdiscrete' and 'rcontinuous', or their multi-step
%            version 'msdiscrete' and 'mscontinuous' with k values ('discrete')
%   - preproc : true to enable preprocessing of inputs ('true')
%   - postproc: [not quite sure yet whether/how it should be used]
%   - ms_k    : multi-step parameter (only for multi-step algorithm)
%
% Output
%   - com: communities (listed for each node)
%   - Q  : stability value of the given partition
%
% Author: Erwan Le Martelot
% Date: 17/02/11

function [com, Q] = gso(adj, ts, model, preproc, postproc, ms_k)
    
    % Set model if not set
    if nargin < 3
        model = 'discrete';
    end
    
    % Enable pre-processing if not specified
    if nargin < 4
        preproc = false;
    end
    
    % Disable post-processing if not specified
    if nargin < 5
        postproc = false;
    end
    
    % Check connections
    single = find(sum(adj)==0);
    if ~isempty(single)
        fprintf('Found single nodes: '); disp(single);
        not_single = setdiff(1:length(adj), single);
        adj = adj(not_single,not_single);
    end
    
    % Pre-process data if required
    if preproc
        % Degree vector
        d = sum(adj,2);
        % Find nodes with a degree of one
        d1 = find(d==1);
        % Get the nodes they should merge with
        if ~isempty(d1)
            %fprintf('Removing nodes: ');
            for i=1:length(d1)
                %fprintf('%d ',d1(i));
                n1(i) = find(adj(d1(i),:)==1);
            end
            n1 = n1';
            %fprintf('\n');
        end
        % Extract adj without the 1 degree nodes
        p_adj = adj(d>1,d>1);
        if sum(sum(p_adj)) == 0
            com = ones(length(adj),1);
            Q = 0;
            return;
        end
    else
        p_adj = adj;
    end
    
    % Call the required algorithm
    switch model
        case {'discrete'}
            [com,Q] = gso_discrete(p_adj, ts);
        case {'rdiscrete'}
            [com,Q] = rgso_discrete(p_adj, ts);
        case {'msdiscrete'}
            if nargin < 6
                [com,Q] = msgso_discrete(p_adj, ts);
            else
                [com,Q] = msgso_discrete(p_adj, ts, ms_k);
            end
        case {'continuous'}
            [com,Q] = gso_continuous(p_adj, ts);
        case {'rcontinuous'}
            [com,Q] = rgso_continuous(p_adj, ts);
        case {'mscontinuous'}
            if nargin < 6
                [com,Q] = msgso_continuous(p_adj, ts);
            else
                [com,Q] = msgso_continuous(p_adj, ts, ms_k);
            end
        otherwise
            error('Unknown model');
    end
    
    % Post-process data if required
    if postproc
        [com, Q] = kernighan_lin(p_adj, com, ts);
    end
    
    % Insert removed nodes back
    if preproc && ~isempty(d1)
        %fprintf('Re-adding nodes: ');
        for i=1:length(d1)
            %fprintf('%d ',d1(i));
            if d1(i) == 1
                com = [0; com];
            else
                com = [com(1:d1(i)-1); 0; com(d1(i):length(com))];
            end
        end
        com(d1) = com(n1);
        %fprintf('\n');
    end
    
    % Add single nodes back with community 0
    if ~isempty(single)
        tmp = [not_single' com; single' zeros(length(single),1)];
        tmp = sortrows(tmp,1);
        com = tmp(:,2);
    end
    
end

