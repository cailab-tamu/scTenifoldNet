function Xhat0=do_td_cp(XM0)
    T0=tensor(XM0);
    [~,U1]=cp_als(T0,5,'printitn',0);
    M2=cp_als(T0,5,'maxiters',100,...
              'init',U1,'printitn',0);
    %Use HOSVD initial guess
    %Use the 'nvecs' option to use the leading mode-n singular vectors as the initial guess.
    %M2 = cp_als(T0,5,'init','nvecs','printitn',0);
    fM0=full(M2);
    Xhat0=fM0.data;
end