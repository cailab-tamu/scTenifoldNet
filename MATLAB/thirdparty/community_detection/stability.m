% Compute the stability of a given partition (Markov chain model)
%
% Input
%   - adj: adjacency matrix
%   - com: community list (per node)
%   - ts : Markov time serie to consider
%
% Output
%   - Q: Stability value
%
% Author: Erwan Le Martelot
% Date: 17/03/11

function [Q] = stability(adj, com, ts)

    % Compute initial indicator matrix
    H = get_indicator_matrix(adj,com);

    % Nodes degree (line) vector
    d = sum(adj,2)';

    % Transition matrix M
    M = diag(d) \ adj;

    % Markov chain transition probabilities (line) vector
    pi = d / sum(d);

    % Markov chain transition probabilities diagonal matrix
    PI = diag(pi);

    % Computation of stability given the time series ts
    Q = 1;
    for s=1:length(ts)
        Rs = H' * ( PI * (M^ts(s)) - pi' * pi ) * H;
        Qs = trace(Rs);
        if (Qs < Q)
            Q = Qs;
        end
    end

end

