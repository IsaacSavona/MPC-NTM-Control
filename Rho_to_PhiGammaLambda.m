function [ Phi, Gamma, Lambda ] =  Rho_to_PhiGammaLambda(Rho1,Rho2,Rho3, A,B,C, kappa,Ts,j_BS,zeta,tau_E,eta_CD,w_dep)
%ABN2PhiGamma Summary of this function goes here
%   Detailed explanation goes here
%   A,B are function, C is a matrix
N  = size(Rho1,2);
% nx = size(B,1);
% nu = size(B,2);
Phi1 = A(Rho1(1),Rho2(1), kappa,Ts,j_BS,zeta,tau_E);
Gamma1 = B(Rho3(1), kappa,Ts,eta_CD,w_dep);
Lambda1 = C;
nx = size(Phi1,1);
nu = size(Gamma1,1);
nu2 = size(Gamma1,2);
nC = size(C,2);

% Phi

Phi = zeros(nx*N,nx);
Phi(1:nx,:) = Phi1;
if N>1
    for j = 2:N
        Phi(j,:) = Phi(j-1,:)*A(Rho1(j),Rho2(j), kappa,Ts,j_BS,zeta,tau_E);
    end
end

% Gamma
Gamma = zeros(nx*N,nu2);
Gamma(1:nx,1:nu2) = Gamma1;
for i = 2:N %rows
    for j = 1:i %columns
        if i ~= j
            %Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = A(Rho1(i-j),Rho2(i-j))*B(Rho3(j));
            Gamma((i-1)*nx+1:nx*i,(j-1)*nu2+1:j*nu2) = A(Rho1(i-j),Rho2(i-j), kappa,Ts,j_BS,zeta,tau_E)*Gamma((i-2)*nx+1:nx*(i-1),(j-1)*nu2+1:j*nu2);
            % Basic idea: new Gamma =
            % newA*previously_Multipled_A's_from_same_column
        else
            Gamma((i-1)*nx+1:nx*i,(j-1)*nu2+1:j*nu2) = B(Rho3(j), kappa,Ts,eta_CD,w_dep);
        end
        %Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = Gamma(j-1,:)*A(Rho1(i-j),Rho2(i-j))*B(Rho3(i));
    end
end
% N=4
% i = 1 --> j = {1}2,3,4 --> i-j = {0}      --> {/}            --> {i-j}
% i = 2 --> j = {1,2}3,4 --> i-j = {1,0}    --> {1,/}          --> {i-j,}
% i = 3 --> j = {1,2,3}4 --> i-j = {2,1,0}  --> {2*1,1,/}      -->
% i = 4 --> j = {1,2,3,4} -> i-j = {3,2,1,0} -> {3*2*1,2*1,1,/} ->

Lambda = zeros(nx*N,nC);
Lambda(1:nx,1:nC) = Lambda1;
for i=2:N
    %(current A)x(previous row) + c
    Lambda((i*nx)-(nx-1):i*nx, 1:nC) = A(Rho1(i),Rho2(i), kappa,Ts,j_BS,zeta,tau_E)*Lambda(((i-1)*nx)-(nx-1):(i-1)*nx, 1:nC) + C ;
end

end


