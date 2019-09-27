load neuron_tensor.mat

for k=1:10
    k
    X(:,:,1,k)=full(Ayng{k});
    X(:,:,2,k)=full(Aold{k});
end
T=tensor(X);
[M1,U1]=cp_als(T,3);
M2=cp_als(T,3,'maxiters',100,'init',U1,'printitn',10);
fM=full(M2);
% a1=fM.data(:,:,1,1);
% a3=fM.data(:,:,3,1);
A0=mean(fM.data(:,:,1,:),4);
A1=mean(fM.data(:,:,2,:),4);


%%
close all

i=400;
%i=2700;
figure;
for k=1:5
subplot(2,6,k)
imagesc(Aold{k}(i:i+49,i:i+49));
title('old')
end
for k=1:5
subplot(2,6,6+k)
imagesc(Ayng{k+1}(i:i+49,i:i+49));
title('yng')
end

subplot(2,6,12)
imagesc(A1(i:i+49,i:i+49))
title('yng denoised')

subplot(2,6,6)
imagesc(A0(i:i+49,i:i+49))
title('old denoised')

