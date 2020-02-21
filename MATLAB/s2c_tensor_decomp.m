% load examplenetworks.mat A0 A1
% clear X
addpath('thirdparty\tensorlab_2016-03-28\');

for k=1:length(Ax0)
    Xx0(:,:,k)=full(Ax0{k});
end
for k=1:length(Ax1)
    Xx1(:,:,k)=full(Ax1{k});
end


Uhat=cpd(Xx0,5);



% A1=tensor_reconst(A,B,C,lam);




A0=mean(A0,3);
A0=A0./max(abs(A0(:)));
A0=round(A0,1);

A1=mean(A1,3);
A1=A1./max(abs(A1(:)));
A1=round(A1,1);
