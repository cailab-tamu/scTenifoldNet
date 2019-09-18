
K = 8; % nearest neighbors for knn
D = 2;  % dimension of manifold
show_corrs = false;

%% synthetic data
N = 600;
[data1,data2] = synthesize_data(N,'fuzzy roll',true);
% data1=sparse(data1);
% data2=sparse(data2);

%% handwritten letter 'b' (uncomment to test)
% data1 = [442 442 442 443 446 439 437 427 419 410 403 397 396 398 404 413 426 439 455 476 499 527 558 589 623 655 683 706 722 730 731 722 707 684 658 629 598 569 540 525 512 508 511 521;
% 523 512 512 524 539 581 637 710 791 888 981 1075 1160 1232 1287 1324 1340 1340 1325 1298 1257 1213 1165 1120 1082 1053 1034 1030 1040 1061 1093 1128 1168 1208 1243 1276 1300 1312 1314 1308 1301 1285 1270 1258];
% data2 = [407 402 402 402 396 393 397 399 404 407 409 412 413 415 414 413 415 415 421 425 433 445 459 480 504 534 569 609 650 691 730 765 792 811 817 811 793 763 722 676 625 572 518 479 442 419 408 408;
% 420 408 408 408 432 472 526 595 677 767 854 941 1017 1087 1135 1168 1181 1181 1171 1151 1118 1084 1041 1004  975  949  934  928  938  954  981 1012 1044 1080 1117 1148 1173 1191 1205 1213 1215 1209 1199 1187 1178 1171 1158 1158];
% data2 = -data2;

%% get stats
[d1,n1] = size(data1);
[d2,n2] = size(data2);

%% unaligned
figure;
plot_correspondences(data1,data2);
title('Unaligned');

%% CCA
if n1 == n2
    [aln1,aln2] = cca(data1,data2,{'b',D});
    figure;
    plot(aln1(1,:),aln1(2,:),'-',aln2(1,:),aln2(2,:),'-');
    title('CCA alignment');
end

%% CTW
if d1 < d2, data1(d1+1:d2,:) = 0; end
if d2 < d1, data2(d2+1:d1,:) = 0; end
Xs = {data1,data2};
[ctwAln, ys] = ctw(Xs, dtw(Xs), [], {'b',D});
data1 = data1(1:d1,:);
data2 = data2(1:d2,:);
aln1 = ys{1};
aln2 = ys{2};

figure;
subplot(1,1+show_corrs,1);
if show_corrs
    plot_correspondences(aln1,aln2,ctwAln.P);
else
    plot_correspondences(aln1,aln2);
end
title('CTW alignment');
if show_corrs
    subplot(1,2,2);
    plot_correspondences(data1,data2,ctwAln.P);
    title('Correspondences in original space');
end

%% MW
disp('>> Linear manifold warping');
[Plin,aln1,aln2] = manifold_warping(data1,data2,'linear',D,K,10);
figure;
subplot(1,1+show_corrs,1);
if show_corrs
    plot_correspondences(aln1,aln2,Plin);
else
    plot_correspondences(aln1,aln2);
end
title('MW linear alignment');
if show_corrs
    subplot(1,2,2);
    plot_correspondences(data1,data2,Plin);
    title('Correspondences in original space');
end

disp('>> Non-Linear manifold warping');
[Pnon,aln1,aln2] = manifold_warping(data1,data2,'nonlinear',D,K,10);
figure;
subplot(1,1+show_corrs,1);
if show_corrs
    plot_correspondences(aln1,aln2,Pnon);
else
    plot_correspondences(aln1,aln2);
end
title('MW nonlinear alignment');
if show_corrs
    subplot(1,2,2);
    plot_correspondences(data1,data2,Pnon);
    title('Correspondences in original space');
end

disp('>> 2-step MW');
figure;
subplot(1,1+show_corrs,1);
[ectwP,aln1,aln2] = manifold_warping(data1,data2,'embed',D,K,10);
title('Two-step MW');
if show_corrs
    plot_correspondences(aln1,aln2,ectwP);
else
    plot_correspondences(aln1,aln2);
end
if show_corrs
    subplot(1,2,2);
    plot_correspondences(data1,data2,ectwP);
    title('Correspondences in original space');
end

%% warping paths
figure;
plot([1, n1], [1, n2], Plin(:,1), Plin(:,2), Pnon(:,1), Pnon(:,2), ...
     ctwAln.P(:,1), ctwAln.P(:,2), ectwP(:,1), ectwP(:,2));
legend('true','MW linear','MW nonlinear','CTW','Two-step MW','Location','Best');
