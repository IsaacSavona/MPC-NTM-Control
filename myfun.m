function OBJ = myfun(u,x0)
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

% Linearized State-Space Formulation
Ts = 0.4;
%A0 = [ , ; , ];
%B0 = [ ; ];
[nx,nu] = size(B0);
N = 5;

% Cost Function Weights
Q = 1*eye(nx);
%R = 0.1*eye(nu);
%[~,P,~] = dlqr(A0,B0,Q,R);

%%% Constructing K and P from LMI %%%
O = sdpvar(nx,nx);
Y = sdpvar(nu,nx);
Con1 = [O, (A*O+B*Y)', O, Y';
    (A*O+B*Y), O, zeros(nx,nx), zeros(nx,nu);
    O, zeros(nx,nx), Q^-1, zeros(nx,nu);
    Y, zeros(nu,nx), zeros(nu,nx), R^-1]>=0;
Con2 = O>=1e-9;
constraints = Con1 + Con2;
diagnostics = optimize(constraints);
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end
P = value(O)^-1;
%K = value(Y)*value(O)^-1;

%------------------------------
x = zeros(nx,N+1);
x(:,1) = x0;
OBJ = 0;
for i = 1:N
    %OBJ = OBJ + x(:,i)'*Q*x(:,i)+u(:,i)'*R*u(:,i);
    OBJ = OBJ + (x(:,i) - r)'*Q*(x(:,i) - r );
    x(:,i+1) = A*x(:,i)+B*u(:,i);
end
OBJ = OBJ + (x(:,N+1) - r)'*P*(x(:,N+1) - r);
end