X = rand(2000,2000);
k = 3;

opts.tol = 1e-8;
opts.maxit = 150;

tic; [U,S,V] = svds(X,k,'L',opts); toc; norm(U*S*V' - X,'fro')
tic; [U,S,V] = lmsvd(X,k,opts); toc; norm(U*S*V' - X,'fro')