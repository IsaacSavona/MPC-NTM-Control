function [ W, L, c, Xi] = getWLcXi(xmax, xmin, umax, umin, Gamma, Phi)
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
      zeros(nx,nu);
%     zeros(nu);
%     -eye(nu); 
%     +eye(nu); %eye dim for u
%     zeros(nu);
%     zeros(nu);
%     zeros(nx*2,nu)]; %make sure the zeros add across given Mi dims

% zetai       = [zeros(nu);
%     zeros(nu);
%     zeros(nu);
%     zeros(nu);
%     -eye(nu);
%     +eye(nu); %eye dim for y
%     zeros(nx*2,nu)]; %make sure the zeros add across given Mi dims

% xii        = [-eye(nu);
%     +eye(nu); %eye dim for du
%     zeros(nu);
%     zeros(nu);
%     zeros(nu);
%     zeros(nu);
%     zeros(nx*2,nu)]; %make sure the zeros add across given Mi dims

bi          = [-umin;
    umax;
    -xmin;
    xmax];

MN = [-eye(nx);eye(nx)];
%zetaN = [-eye(nu);eye(nu);zeros(nx,nu);zeros(nx,nu)];
bN = [-xmin;xmax];

%zetaN = [-ymin;xmax];
% uncomment previous 2 lines to not use the terminal set constraints

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

% %%% cal H %%%
% Hcal = [zetai;repmat(0*zetai,N-1,1);0*zetaN];
% 
% %%% Zeta %%%
% Zeta = zetaN;
% for i = 2:N
%     Zeta = blkdiag(zetai,Zeta);
% end
% Zeta = [zeros(size(zetai,1),size(Zeta,2));Zeta];
% 
% %%% Xi %%%
% Xi = xii;
% for i = 2:N
%     Xi = blkdiag(Xi,xii);
% end
% Xi = [Xi;zeros(size(zetaN,1),size(Xi,2))];

c = bN;
for i = 1:N
    c = [bi;c];
end

%L = Mcal*Gamma+Ecal;
%size(Mcal*Gamma)
%size(Zeta)
%size(Lambda*Gamma)
%size(Zeta*Lambda*Gamma)
%size(Ecal)
%size(Xi)
%size(Xi*Is)
%Y = Mcal*Gamma + Zeta*Lambda*Gamma;
L =Mcal*Gamma + Ecal
W = -Dcal-Mcal*Phi;
W = -Dcal - Hcal*C - Mcal*Phi - Zeta*Lambda*Phi;
%W = W;
%Xi = -Xi;

end