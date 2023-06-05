function B=B(rho3, wdep, kappa, Ts, etaCD)
    global ;
    B=[(kappa*Ts*etaCD/wdep)*rho3 0];
end