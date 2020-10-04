function Xhat0=do_td_tucker(XM0)
    T0=tensor(XM0);
    M2=tucker_als(T0,5,'printitn',0);
    fM0=full(M2);
    Xhat0=fM0.data;
end