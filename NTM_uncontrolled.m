%% Reset
clear all

%% Physics parameters
w_dep = 0.024; % [m] deposition width
tau_r = 293;   % [s] resistive time scale
rs = 1.55;     % [m] radial location of island
a = 2.0;       % [m] minor radius 
w_marg = 0.02; % [m] critical width
eta_CD = 0.9;  % [#] current drive efficiency
tau_E0 = 3.7;  % [s] energy confinement time without islands
tau_E = tau_E0;% [s] currently NOT EXACT FORMULA
mu0 = 4e-7*pi; % [?] 
Lq = 0.87;     % [m] q-gradient length scale
B_pol = 0.97;  % [T] poloidal field
m = 2;         % [#] poloidal mode number
Cw = 1;        % [#] UNKNOWN
tau_A0 = 3e-6; % [s] Alfvèn time
tau_w = 0.188; % [s] resistive wall time
kappa = 16*mu0*Lq*rs^2/(0.82*tau_r*B_pol*pi);
zeta = m*Cw*tau_A0^2*tau_w;


%% Quasi-LPV MPC Model %%

%%% Model
N = 3;    % prediction horizon
Ts = 0.1; % sampling time
nx = 2;   % dimensions of state vector
nu = 1;   % dimensions of input vector
x0 = [0.01;1000]; % initial state (low width and high frequency does not need control)

%%% Compact Formulation
Rho1 = repmat(rho1(x0(:),w_marg),1,N); % initial Rho by using current rho(k) at every predicted step
Rho2 = repmat(rho2(x0(:)),1,N);

%% Quasi-LPV MPC Simulation %%

%%% initialize variables
k_sim = 20; % number of simulation time steps
xk = [x0 zeros(nx,k_sim)]; % states [w(k),ω(k)] at every time step k=0 ... k=k_sim (size=(nx) x (k_sim+1))

%%% Controller Iterations
for k = 1:k_sim              % simulation loop over time samples
    xk(:,k+1) = A(rho1(xk(:,k),w_marg), rho2(xk(:,k)), kappa, rs, tau_r, Ts, zeta, a, tau_E)*xk(:,k); % evolve state one time step
end


%% Plot 4.1.1

%%% Plots the states and input trajectories over time%%%
figure
subplot(1,1,1)
xkPlot = xk';
hold on
yyaxis left
plot(0:k_sim,xk(1,:))
ylabel('$\mathrm{w}$ [m]','Interpreter','latex')
yyaxis right
plot(0:k_sim,xk(2,:))
hold off
xlabel('$k$','Interpreter','latex')
ylabel('$\omega$ [Hz]','Interpreter','latex')
% title('Constrained quasi-LPV MPC State Trajectory')
%axis([0 k_sim -5 10])




%% TEMPORARY functions
function rho1 = rho1(x,w_marg)
    rho1 = 1/(x(1)^2 + w_marg^2);
end

function rho2 = rho2(x)
    rho2 = x(1)^2/x(2);
end

function A = A(rho1, rho2, kappa, rs, tau_r, Ts, zeta, a, tau_E)
    A = [((4/3)*(kappa*rs/(0.82*tau_r))*Ts*rho1+1)  0; 
        ((rho2*Ts)/(zeta*a^3)) (1-Ts/tau_E)];
end