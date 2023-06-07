function rho3 = rho3(x, w_dep)
    wstar = x(1)/w_dep;
    rho3 = (0.25 + 0.24*wstar)/(1 + 1.5*wstar + 0.43*wstar^2 + 0.64*wstar^3);
end
