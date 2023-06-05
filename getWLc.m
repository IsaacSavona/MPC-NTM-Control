function [W, L, c] = getWLc(xmax, xmin, umax, umin, Gamma, Phi)
%GETWLC Summary of this function goes here
%   Detailed explanation goes here

nu = size(umin,1);
nx = size(xmin,1);
N = size(Phi,1)/nx;

Mi          = [zeros(nu,nx);
    zeros(nu,nx);
    -eye(nx) %eye dim for x
    +eye(nx)];

Ei          = [-1;
      1;
      zeros(nx,nu);
      zeros(nx,nu)];


bi          = [-umin;
    umax;
    -xmin;
    xmax];

MN = [-eye(nx);eye(nx)];
bN = [-xmin;xmax];


%%% cal D %%%
Dcal = [Mi;repmat(0*Mi,N-1,1);0*MN];

%%% cal M %%%
Mcal = MN;
for i = 2:N
    Mcal = blkdiag(Mi,Mcal);
end
Mcal = [zeros(size(Mi,1),size(Mcal,2));Mcal];

%%% cal E %%%
Ecal = Ei;
for i = 2:N
    Ecal = blkdiag(Ecal,Ei);
end
Ecal = [Ecal;zeros(size(MN,1),size(Ecal,2))];

c = bN;
for i = 1:N
    c = [bi;c];
end

L = Mcal*Gamma + Ecal;
W = -Dcal - Mcal*Phi;


end