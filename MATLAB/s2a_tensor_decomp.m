% load examplenetworks.mat A0 A1
% clear X
addpath('thirdparty\tensor_toolbox-v3.1\');
for k=1:length(Ax0)
    Xx0(:,:,k)=full(Ax0{k});
end
for k=1:length(Ax1)
    Xx1(:,:,k)=full(Ax1{k});
end

T0=tensor(Xx0);
[M1,U1]=cp_als(T0,3);
M2=cp_als(T0,3,'maxiters',100,'init',U1,'printitn',10);
fM0=full(M2);

T1=tensor(Xx1);
[M1,U1]=cp_als(T1,3);
M2=cp_als(T1,3,'maxiters',100,'init',U1,'printitn',10);
fM1=full(M2);

A0=mean(fM0.data,3);
A1=mean(fM1.data,3);

A0=A0./max(abs(A0(:)));
A0=round(A0,1);
A1=A1./max(abs(A1(:)));
A1=round(A1,1);

%A0=A0-diag(diag(A0));
%A1=A1-diag(diag(A1));

%A0=A0.*(abs(A0)>quantile(abs(A0(:)),0.95));
%A1=A1.*(abs(A1)>quantile(abs(A1(:)),0.95));

% writematrix(genelist,'genelist.txt')
% writematrix(A0,'A0');
% writematrix(A1,'A1');

%{
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
%}