load neuron_tensor.mat
for k=1:10
    X(:,:,1,k)=full(Ayng{k});
    X(:,:,2,k)=full(Aold{k});
end
T=tensor(X);
[M1,U1]=cp_als(T,3);
M2=cp_als(T,3,'maxiters',100,'init',U1,'printitn',10);
fM=full(M2);
A0=mean(fM.data(:,:,1,:),4);
A1=mean(fM.data(:,:,2,:),4);

%%
close all

% i=4100;
i=400;
ix=i:i+50;

%i=2700;
figure;
for k=1:5
    subplot(2,6,k)
    imagesc(Aold{k}(ix,ix));
    title('old')
    subplot(2,6,6+k)
    imagesc(Ayng{k+1}(ix,ix));
    title('yng')
end
subplot(2,6,12)
imagesc(A1(ix,ix))
title('yng denoised')
subplot(2,6,6)
imagesc(A0(ix,ix))
title('old denoised')

% figure;
% for k=1:5
%     subplot(2,6,k)
%     a=fM.data(:,:,1,k);
%     imagesc(a(ix,ix));
%     title('old')
%     subplot(2,6,6+k)
%     a=fM.data(:,:,2,k);
%     imagesc(a(ix,ix));
%     title('yng')
% end

% clearvars -except genelist A0 A1
