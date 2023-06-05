function Gamma, Phi = GammaPhi(rho1, rho2, rho3)
    global w_dep kappa taur Ts zeta rs;
   
    N = length(rho1)
    Gamma = zeros(N,N;
    Phi = zeros(N);
    A = A(rho1, rho2);
    B = B(rho3);

    for i=1:N %index for rows
        Phi(i) = prod(A(i:N);
        for ii=1:N %index for columns
            Gamma(i,ii) = prod(A(1:i-1),rho2(i-1)) * prod(B(0:ii));
        end
    end
end