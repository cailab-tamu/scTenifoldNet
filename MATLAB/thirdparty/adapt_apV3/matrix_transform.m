function [M,dim] = matrix_transform(M,ngiven,nmax)

Ms = -realmax*ones(ngiven,ngiven);
for i=1:ngiven
    Ms(i,i) = 0;
end

% transform to full matrix
if 1 % elements from partial matrix
    ns = size(M,1);
    for i = 1:ns
        ni = M(i,1);
        nj = M(i,2);
        if (ni < ngiven || ni == nmax) && (nj < ngiven || nj == nmax)
            if ni == nmax
                ni = ngiven;
            end
            if nj == nmax
                nj = ngiven;
            end
            Ms(ni,nj) = M(i,3);
        end
    end
else % elements from full matrix
    m = ngiven-1;
    ns = 1:m;
    for i = 1:3
        ns = [ns i*(nmax-1)+(1:m)];
    end
    for i = ns
        ni = M(i,1);
        nj = M(i,2);
        if ni <= ngiven && nj <= ngiven
            Ms(ni,nj) = M(i,3);
        end
    end
end

% transform full matrix to vectors
nap = ngiven*ngiven-ngiven;
M = zeros(nap,3);
j=1;

for i = 1:ngiven
   for k = [1:i-1,i+1:ngiven]
     M(j,1) = i;
     M(j,2) = k; 
     M(j,3) = Ms(i,k);
     j = j+1;
   end;
 end;
dim = length(M);
