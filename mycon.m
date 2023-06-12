function [CON,CONe] = mycon(u,x0,M_N,b_N)
%------------------------------
% Physics Parameters
j_BS = 73e3;   % [A/m^2] bootstrap current density
w_dep = 0.024; % [m] deposition width
w_marg = 0.02; % [m] marginal width
w_sat = 0.32;  % [m] saturation width
tau_r = 293;   % [s] resistive time scale
rs = 1.55;     % [m] radial location of island
a = 2.0;       % [m] minor radius 
eta_CD = 0.9;  % [#] current drive efficiency
tau_E0 = 3.7;  % [s] energy confinement time without islands
tau_E = tau_E0;% [s] currently NOT EXACT FORMULA
mu0 = 4e-7*pi; % [N/A^2] vacuum permaebility
Lq = 0.87;     % [m] q-gradient length scale
B_pol = 0.97;  % [T] poloidal field
m = 2;         % [#] poloidal mode number
Cw = 1;        % [#] UNKNOWN!!
tau_A0 = 3e-6; % [s] Alfvèn time
tau_w = 0.188; % [s] resistive wall time
omega0 = 2*pi*420; % [rad/s] equilibrium frequency

%%% Non-Linear State-Space
N = 20; % Make sure this is consistent with NMPC_Sim
Ts = 0.4;
%A = [ , ; , ]; Nonlinear A
%B = [ ; ]; Nonlinear B
[nx,nu] = size(B);

%%% Constraints
min_width = 0.06; % [m] RANDOM VALUE based on 2.4cm deposition width + error in dep. pos.
max_width = 0.15; % [m]
min_freq = 100*2*pi;  % [rad/s]
max_freq = 5000*2*pi; % [rad/s]
xmin = [min_width;min_freq]; % minimum on state vector
xmax = [max_width;max_freq]; % maximum on state vector

min_power = 0;    % [W]
max_power = 2e6;  % [W]
umin = min_power; % minimum on input vector
umax = max_power; % maximum on input vector
%------------------------------
CON = [];
CONe = [];
x = zeros(nx,N+1);
x(:,1) = x0;
for i = 1:N
    x(:,i+1) = A*x(:,i)+B*u(:,i);
    CON = [CON; -x(:,i+1)+xmin; x(:,i+1)-xmax; -u(:,i)+umin; u(:,i)-umax];
end
CON = [CON; M_N*x(:,N+1)-b_N];
end