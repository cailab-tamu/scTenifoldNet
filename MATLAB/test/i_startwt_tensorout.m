A0=csvread('tensorOut_X.csv',1,1);
A1=csvread('tensorOut_Y.csv',1,1);
run('../s3_manifold_algn')
[~,ind]=sort(vecnorm(aln0-aln1,2,2),'descend');
ind(1:5)