
clc
clear

import casadi.*
%%
global nx nu nd
global lbx ubx dx0 lbu ubu u0

par.tf = 0.1;
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
x0 = [0.04;5000*2*pi]; % initial state (low width and high frequency does not need control)

%%

sys.N = 20;
sys.discrete = 0;

%sys.d = 0; % model disturbance 
%plant.d= 2; % plant disturbance
%plant.d= 0; % plant disturbance

sys.P = eye(nx); % terminal cost
%sys.X_N = [0;0];    % terminal state
r = [0.10*1e2; 1000*2*pi*1e-3]; % reference state
%sys.X_N = r;    % terminal state
ADP = BuildADP_N(sys);

%%

%x =[pi;0]; % simulation initial condition
x = [0.04*1e2;5000*2*pi*1e-3]; % simulation initial condition
sim.V = 0;
Primal = ADP.x0;
for i = 1:120/par.tf

    tic;
    sol = ADP.solver('x0',Primal,'p',vertcat(x,r),...
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
figure(12)
subplot(2,2,1)
hold all
plot((1:1:120/par.tf)*par.tf,sim.x(:,1))
subplot(2,2,3)
hold all
stairs((1:1:120/par.tf)*par.tf,sim.u)
subplot(2,2,[2,4])
hold all
plot(sim.sol_t)

