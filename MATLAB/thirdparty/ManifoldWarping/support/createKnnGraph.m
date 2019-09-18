function W = createKnnGraph(N1)
	
n1 = size(N1,1);
W1 = sparse(n1,n1);
%TODO: vectorize?
for i=1:n1
    W1(i,N1(i,:)) = 1;
end
W = (W1+W1')/2;

end
