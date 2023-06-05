function [ Phi, Gamma ] =  Rho_to_PhiGamma(Rho1,Rho2,Rho3, A, B)
%ABN2PhiGamma Summary of this function goes here
%   Detailed explanation goes here

N  = size(Rho1,1); 
% nx = size(B,1);
% nu = size(B,2);
Phi1 = A(Rho1(1),Rho2(1));
Gamma1 = B(Rho3(1));
nx = size(Phi1,1);
nu = size(Gamma1,1);

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
Gamma(1:nx,1:nu) = B(Rho3(1));
for i = 2:N
    for j = 1:i
        if i == j
            %Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = A(Rho1(i-j),Rho2(i-j))*B(Rho3(j));
            Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = A(Rho1(i-j),Rho2(i-j))*Gamma((i-2)*nx+1:nx*(i-1),(j-1)*nu+1:j*nu);
            % Basic idea: new Gamma = newA*previously_Multipled_A's
        else
            Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = B(Rho3(j));
        end
        %Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = Gamma(j-1,:)*A(Rho1(i-j),Rho2(i-j))*B(Rho3(i));
    end
end
% N=4
% i = 1 --> j = {1}2,3,4 --> i-j = {0}      --> {/}            --> {i-j}
% i = 2 --> j = {1,2}3,4 --> i-j = {1,0}    --> {1,/}          --> {i-j,}
% i = 3 --> j = {1,2,3}4 --> i-j = {2,1,0}  --> {2*1,1,/}      -->
% i = 4 --> j = {1,2,3,4} -> i-j = {3,2,1,0} -> {3*2*1,2*1,1,/} ->
end

