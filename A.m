function A=A(rho1, rho2, kappa, Ts, j_BS, zeta, tau_E)
    A=[(1+(4/3)*(kappa*Ts*j_BS*rho1))  0; (-(rho2*Ts)/(zeta)) (1-Ts/tau_E)];
end
