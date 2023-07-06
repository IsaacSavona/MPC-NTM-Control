
clc
clear

import casadi.*
%%
global nx nu nd
global lbx ubx dx0 lbu ubu u0

par.tf = 0.1;
par.T = 20;
[sys,F_integrator] = NTM(par.tf);

nx = numel(sys.x);
nu = numel(sys.u);
nd = 0;

%lbu = -20*ones(nu,1);
%ubu = 20*ones(nu,1);
lbu = 1e-3.*ones(nu,1);  %[MW] min power
ubu = 2*ones(nu,1); %[MW] max power
u0  = 1;

% lbx = [0-pi,-5]';
% ubx = [2*pi-pi,5]';
min_width = 0.00; % [m] RANDOM VALUE based on 2.4cm deposition width + error in dep. pos.
max_width = 0.15; % [m]
min_freq = 100*2*pi;  % [rad/s]
max_freq = 5000*2*pi; % [rad/s]
lbx = [min_width*1e2,min_freq*1e-3]';
ubx = [max_width*1e2,max_freq*1e-3]';
%dx0 = [0;0];

%%

sys.N = 20;
sys.discrete = 0;

%sys.d = 0; % model disturbance 
%plant.d= 2; % plant disturbance
%plant.d= 0; % plant disturbance

sys.P = eye(nx); % terminal cost
%sys.X_N = [0;0];    % terminal state
%r = [0.002*1e2; 1000*2*pi*1e-3]; % reference state
k_sim = par.T/par.tf;
r = ones(nx,k_sim).*[1e2;1e-3];
r(:,1:round(k_sim*(1/5))) = r(:,1:round(k_sim*(1/5))).*[0.04; 5000*2*pi];
r(:,round(k_sim*(1/5))+1:round(k_sim*(1/5+1/6))) = r(:,round(k_sim*(1/5))+1:round(k_sim*(1/5+1/6))).*[0.005; 5000*2*pi];
r(:,round(k_sim*(1/5+1/6))+1:round(k_sim*(1/5+1/6+1/4))) = r(:,round(k_sim*(1/5+1/6))+1:round(k_sim*(1/5+1/6+1/4))).*[0.03; 5000*2*pi];
r(:,round(k_sim*(1/5+1/6+1/4))+1:round(k_sim*(1/5+1/6+1/2))) = r(:,round(k_sim*(1/5+1/6+1/4))+1:round(k_sim*(1/5+1/6+1/2))).*[0.02; 5000*2*pi];
r(:,round(k_sim*(1/5+1/6+1/2))+1:k_sim) = r(:,round(k_sim*(1/5+1/6+1/2))+1:k_sim).*[0.10; 5000*2*pi];

%sys.X_N = r;    % terminal state
ADP = BuildADP_N(sys);

%%

%x =[pi;0]; % simulation initial condition
x = [0.04*1e2;5000*2*pi*1e-3]; % simulation initial condition
sim.V = 0;
Primal = ADP.x0;
for i = 1:par.T/par.tf

    tic;
    sol = ADP.solver('x0',Primal,'p',vertcat(x,r(:,1)),...
        'lbx',ADP.lbx,'ubx',ADP.ubx,'lbg',ADP.lbg,'ubg',ADP.ubg);
    sim.sol_t(i) = toc;
    
    flag = ADP.solver.stats();
    assert(flag.success, ['ADP solver failed: ' flag.return_status])
    
    Primal = full(sol.x); 
    u = Primal(1); % optimal input
    sim.x(i,:) = x.*[1e-2;1e3];% Scaling back to the orignial units for plotting
    sim.u(i) = u*1e6; % Scaling back to the orignial units for plotting
    sim.t(i) = i*par.tf;
    
    % --------- Nonlinear plant simulation -------------
    Fk = F_integrator('x0',x,'p',u);
    x =  full(Fk.xf);
end
%%
figure(9)
subplot(2,2,1)
hold all
plot((1:1:par.T/par.tf)*par.tf,sim.x(:,1)*1e2)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')

subplot(2,2,3)
hold all
stairs((1:1:par.T/par.tf)*par.tf,sim.u/1e6)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$P_{ECCD}$ [MW]','Interpreter','latex')

subplot(2,2,[2,4])
semilogy(sim.sol_t*1e3)
hold all
xlabel('sample $k$ [\#]','Interpreter','latex')
ylabel('computation time $\tau$ [ms]','Interpreter','latex')



