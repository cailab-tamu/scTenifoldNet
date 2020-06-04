load X.txt
[A]=sc_pcnet(X);
A=A./max(abs(A(:)));
A(1:5,1:5)
