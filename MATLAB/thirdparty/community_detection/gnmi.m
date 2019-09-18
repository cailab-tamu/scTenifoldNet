% Compute the generalised normalised mutual information between two
% sets of overlapping communities as defined in (Lancichinetti, Fortunato,
% Kertész; Detecting the overlapping and hierarchical community structure
% in complex networks, New Journal of Physics, 2009)
%
% Input:
%   - C1: first set of communities
%   - C2: second set of communities
%   - N : number of nodes in the network
%
% Output:
%   - NXY: normalised mutual information
%
% Author: Erwan Le Martelot
% Date: 21/11/10

function [NXY] = gnmi(C1, C2, N)

    % H(X)
    for c=1:length(C1)
        p = length(C1{c})/N;
        if (p == 0) || (p == 1)
            HX(c) = 0;
        else
            HX(c) = -p*log(p) - (1-p)*log(1-p);
        end
    end
    
    % H(Y)
    for c=1:length(C2)
        p = length(C2{c})/N;
        if (p == 0) || (p == 1)
            HY(c) = 0;
        else
            HY(c) = -p*log(p) - (1-p)*log(1-p);
        end
    end

    % H(Xk|Yl)
    nbtestc1 = zeros(length(C1),1);
    minHXY = Inf(length(C1),1);
    nbtestc2 = zeros(length(C2),1);
    minHYX = Inf(length(C2),1);
    for c1=1:length(C1)
        lc1 = length(C1{c1});
        for c2=1:length(C2)
            lc2 = length(C2{c2});
            %l12 = length(intersect(C1{c1},C2{c2}));
            l12 = intersect_size(C1{c1},C2{c2});
            p11 = l12/N;
            p10 = (lc1 - l12)/N;
            p01 = (lc2 - l12)/N;
            p00 = (N - (lc1 + lc2 - l12))/N;
            if p11 > 0
                h11 = - p11 * log(p11);
            else
                h11 = 0;
            end
            if p10 > 0
                h10 = - p10 * log(p10);
            else
                h10 = 0;
            end
            if p01 > 0
                h01 = - p01 * log(p01);
            else
                h01 = 0;
            end
            if p00 > 0
                h00 = - p00 * log(p00);
            else
                h00 = 0;
            end
            HXY = (h11 + h10 + h01 + h00) - HY(c2);
            HYX = (h11 + h10 + h01 + h00) - HX(c1);
            if (h11 + h00) > (h01 + h10)
                nbtestc1(c1) = nbtestc1(c1) + 1;
                if HXY < minHXY(c1)
                    minHXY(c1) = HXY;
                end
                nbtestc2(c2) = nbtestc2(c2) + 1;
                if HYX < minHYX(c2)
                    minHYX(c2) = HYX;
                end
            end
        end
    end
    
    % H(Xk|Y) norm summed
    HXYn = 0;
    for c=1:length(C1)
        if (nbtestc1(c)>0) && (HX(c) > 0)
            HXYn = HXYn + minHXY(c) / HX(c);
        else
            HXYn = HXYn + 1;
        end
    end
    HXYn = HXYn/length(C1);
    
    % H(Yk|X) norm summed
    HYXn = 0;
    for c=1:length(C2)
        if (nbtestc2(c) > 0) && (HY(c) > 0)
            HYXn = HYXn + minHYX(c) / HY(c);
        else
            HYXn = HYXn + 1;
        end
    end
    HYXn = HYXn/length(C2);
    
    % N(X|Y)
    NXY = 1 - (HXYn + HYXn)/2;

end

% Size of intersection (sets need to be sorted)
function [counter] = intersect_size(c1, c2)
    counter = 0;
	i = 1;
    j = 1;
	while (i<=length(c1)) && (j<=length(c2))
		while (i<=length(c1)) && (c1(i) < c2(j)) i = i+1; end
		if i > length(c1) break; end
		while (j<=length(c2)) && (c2(j) < c1(i)) j = j+1; end
		if j > length(c2) break; end
		if c1(i) == c2(j)
            j = j+1;
            counter = counter+1;
        end
        i = i+1;
    end
end
