% Greedy stabilisation optimisation with overlapping communities
% given a model with possibility of input pre-processing and output
% post-processing
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
%   - ecom: edge communities (listed for edge i,j)
%   - ncom: node communities (listed for each node)
%   - Q   : stability value of the given partition
%
% Author: Erwan Le Martelot
% Date: 18/05/11

function [ecom,ncom,Q] = ogso(adj, ts, model, preproc, postproc, ms_k)
    [lg_adj, lg_nodes] = line_graph(adj);
    if nargin < 3
        [lg_com, Q] = gso(lg_adj, ts);
    elseif nargin < 4
        [lg_com, Q] = gso(lg_adj, ts, model);
    elseif nargin < 6
        [lg_com, Q] = gso(lg_adj, ts, model, preproc, postproc);
    else
        [lg_com, Q] = gso(lg_adj, ts, model, preproc, postproc, ms_k);
    end
    ncom = lgcom2gcom(adj, lg_nodes, lg_com);
    ecom = [lg_nodes lg_com];
end
