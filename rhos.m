w_dep = 1;
N = 2;
kappa = 3;
taur = 1;
Ts = 1;
zeta = 1;
rs = 1;
a = 2.8;
w_marg = 2e-2;



Hello = Gamma([1 1], [2 2], [3 3]);



function rho1 = rho1(x)
    rho1 = 1/(x(1)^2 + w_marg^2);
end

function rho2 = rho2(x)
    rho2 = x(1)^2/x(2)
end

function rho3 = rho3(x)
    wstar = x(1)/w_dep
    rho3 = (0.25 + 0.24*wstar)/(1 + 1.5*wstar + 0.43*wstar^2 + 0.64*wstar^3)
end


function [Gamma, Phi] = Rho_to_PhiGamma(rho1, rho2, rho3)
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

function [W, L, c] = getWLc(N, ulow, uhigh, xlow, xhigh)
    
    %calculate L
    Mu = zeros(N,N);
    for i=2:N-1
        Mu(i,i) = [0 0; -eye(2) eye(2)];
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




%%%%%%%%%%%%%%%%%%%%%%%%% from Cas
function B=B(rho3)
    global N w_dep kappa taur Ts zeta rs;
    B=[(kappa*Ts*etaCD/wdep)*rho3 0];
end

function A=A(rho1, rho2)
    global N w_dep kappa taur Ts zeta rs a;
    A=[((4/3)*(kappa*rs/(0.82*taur))*Ts*rho1+1)  0; ((rho2*Ts)/(zeta*a^3)) (1-Ts/TE)];
end
