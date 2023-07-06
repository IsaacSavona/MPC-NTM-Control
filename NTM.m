function [sys,F] = NTM(Ts)
import casadi.*
% Physics parameters
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
Cw = 1e3;      % [#] UNKNOWN!!
tau_A0 = 3e-6; % [s] Alfvèn time
tau_w = 0.188; % [s] resistive wall time
omega0 = 2*pi*420; % [rad/s] equilibrium frequency

kappa = 16*mu0*Lq*rs^2/(0.82*tau_r*B_pol*pi);
zeta = m*Cw*tau_A0^2*tau_w*a^3;

x = MX.sym('x',2);
u = MX.sym('u',1);
%d = MX.sym('d',1);

dx1 = ((4/3)*(kappa*j_BS*(1/((x(1)*0.01)^2+w_marg^2))))*(x(1)*0.01)...
- (kappa*eta_CD/w_dep)*((0.25 + 0.24*(x(1)*0.01)/w_dep)/(1 + 1.5*(x(1)*0.01)/w_dep + 0.43*((x(1)*0.01)/w_dep)^2 +0.64*((x(1)*0.01)/w_dep)^3))*u(1).*1e6 -...
4/3*(kappa*j_BS*w_sat)/(w_sat^2+w_marg^2); %3.1
dx2 = -(1/zeta)*((x(1)*0.01)^2/(x(2)*1e3))*(x(1)*0.01)...
    - (1/tau_E)*(x(2)*1e3)...
    + omega0/tau_E0; %3.2


dx = vertcat(dx1,dx2);
%rw = 0.001; % [m] reference width

%L = ((x(1)*0.01)).^2 + ((x(2)*1e3)).^2 + 0.01*u.^2;
%L = ((x(1)*0.01)-rw).^2; % punish distance from reference
L = 0;

% create CVODES integrator
ode = struct('x',x,'p',vertcat(u),'ode',dx,'quad',L);
%ode = struct('x',x,'ode',dx,'quad',L);
opts = struct('tf',Ts);
F = integrator('F','cvodes',ode,opts);

f = Function('f',{x,u},{dx});
%f = Function('f',{x,u},{dx,L},{'x'},{'xdot','qj'});

% Runge-Kutta 4 integration
k1 = f(x,u);
k2 = f(x+Ts/2*k1,u);
k3 = f(x+Ts/2*k2,u);
k4 = f(x+Ts*k3,u);
x_next = x+Ts/6*(k1+2*k2+2*k3+k4); 

f_discrete = Function('f_discrete',{x,u},{x_next}); % Discretized model
sys = struct('x',x,'u',u,'dx',dx,'f',f,'F',f_discrete);

sys.D = [1;-1];
%sys.D = [1;1];
sys.bu=[2;0];
