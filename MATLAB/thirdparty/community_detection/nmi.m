% Compute the normalised mutual information between two lists
% Depending on the sizes of the lists the implementation will use a full
% or a sparse representation to optimise speed and memory allocation.
%
% Input:
%   - l1: first list of values
%   - l2: second list of values
%
% Output:
%   - v: normalised mutual information between the two lists
%
% Author: Erwan Le Martelot
% Date: 17/11/10

function [ v ] = nmi( l1, l2 )

    if length(l1) ~= length(l2)
        error('The two lists have different lengths');
    end

    % Number of nodes
    n = length(l1);

    % Number of values in each list
    nc1 = length(unique(l1));
    nc2 = length(unique(l2));

    % Terms initialisation
    term1 = 0;
    term2 = 0;
    term3 = 0;
    
    use_full_matrix = (nc1*nc2 < 1000000);
    % Computing terms using full matrix representation
    if use_full_matrix
        
        % Build the confusion matrix
        c = zeros(nc1,nc2);
        for i=1:n
            c(l1(i),l2(i)) = c(l1(i),l2(i)) + 1;
        end
        sumci = sum(c,2);
        sumcj = sum(c,1);

        % Terms computing
        for i=1:nc1
            for j=1:nc2   
                if c(i,j) > 0
                    term1 = term1 + ( c(i,j) * log( (c(i,j)*n) / (sumci(i)*sumcj(j)) ) );
                end
            end
            term2 = term2 + ( sumci(i) * log(sumci(i)/n) );
        end
        for j=1:nc2
            term3 = term3 + ( sumcj(j) * log(sumcj(j)/n) );
        end
        
    % Sparse representation
    else
        
        % Build the confusion matrix
        c = sparse(nc1,nc2);
        sumci = zeros(nc1,1);
        sumcj = zeros(nc2,1);
        for i=1:n
            c(l1(i),l2(i)) = c(l1(i),l2(i)) + 1;
            sumci(l1(i)) = sumci(l1(i)) + 1;
            sumcj(l2(i)) = sumcj(l2(i)) + 1;
        end

        % Terms computing
        for i=1:nc1
            cols = find(c(i,:));
            for k=1:length(cols)
                j = cols(k);
                term1 = term1 + ( c(i,j) * log( (c(i,j)*n) / (sumci(i)*sumcj(j)) ) );
            end
            term2 = term2 + ( sumci(i) * log(sumci(i)/n) );
        end
        for j=1:nc2
            term3 = term3 + ( sumcj(j) * log(sumcj(j)/n) );
        end
        
    end
     
    % Result
    v = (-2 * term1) / (term2 + term3);

end

