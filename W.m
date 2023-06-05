function W, L, c = WLc(N, ulow, uhigh, xlow, xhigh)
    
    %calculate L
    Mu = zeros(N,N);
    for i=2:N-1
        Mu(i,i) = [0 0 -eye(2) eye(2)];
    end
    Mu(N,N) = [-eye(2) eye(2)];

    Epsilon = zeros(N,N);
    for i=1:N-1
        Epsilon(i,i) = [-1 1 0 0];
    end   
    
    L = Mu*Gamma(rho1,rho2,rho3) + Epsilon;
    %
    
    %calculate W
    D = zeros(N);
    D(1) = [0 0 -eye(2) eye(2)];
    W = -D - Mu*Phi(inputs);
    %

    %calculate c
    c = zeros(N);
    for i = 1:N-1
        c(i) = [-ulow uhigh -xlow xhigh];
    end
    c(N) = [-xlow xhigh];
    %
end