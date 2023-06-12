function B = B(rho3, kappa, Ts, eta_CD, w_dep)
    B = [-(kappa*Ts*eta_CD/w_dep)*rho3; 0]*1e6;
end
