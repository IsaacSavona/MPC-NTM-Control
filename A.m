function A=A(rho1, rho2, kappa, taur, Ts, zeta, rs, a, TE)
    A=[((4/3)*(kappa*rs/(0.82*taur))*Ts*rho1+1)  0; ((rho2*Ts)/(zeta*a^3)) (1-Ts/TE)];
end
