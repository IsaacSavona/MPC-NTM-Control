function [ Phi, Gamma ] =  Rho_to_PhiGamma(Rho1,Rho2,Rho3, A, B)
%ABN2PhiGamma Summary of this function goes here
%   Detailed explanation goes here

N  = size(Rho1)
% nx = size(B,1);
% nu = size(B,2);
Phi1 = A(Rho1(1),Rho2(1));
Gamma1 = B(Rho3(1))
nx = size(Phi1,1)
nu = size(Gamma1,1)

% Phi

Phi = zeros(nx*N,nx);
Phi(1,:) = Phi1;
if N>1
    for j = 2:N
        %Phi     = [Phi;A^j];
        Phi(j,:) = Phi(j-1,:)*A(Rho1(j),Rho2(j));
    end
end

% Gamma
Gamma = zeros(nx*N,nu*N);

for i = 1:N
    for j = 1:i
        Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = A^(i-j)*B;
    end
end
end

