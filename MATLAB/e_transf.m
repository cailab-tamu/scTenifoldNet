function a=e_transf(a)
    a=a./max(abs(a(:)));
    a=a.*(abs(a)>quantile(abs(a(:)),0.95));
end