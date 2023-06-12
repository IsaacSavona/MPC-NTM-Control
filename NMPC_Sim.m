%% Reset
clear all
close all

%% Physics parameters
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
kappa = 16*mu0*Lq*rs^2/(0.82*tau_r*B_pol*pi);
zeta = m*Cw*tau_A0^2*tau_w*a^3;
Ts = 0.4;

Const = [-4/3*(kappa*Ts*j_BS*w_sat)/(w_sat^2+w_marg^2); Ts*omega0/tau_E0];

%% Non-Linear Discrete Model
% A = [(1+(4/3)*(kappa*Ts*j_BS*(rho1(x(1),w_marg))))+(C(1,:)/x(1))  0; (-(rho2(x)*Ts)/(zeta))+(C(1,:)/x(1)) (1-Ts/tau_E)];
% B = [-(kappa*Ts*eta_CD/w_dep)*rho3(x, w_dep); 0];

%% Linearized Discrete Model around w_crit
N = 20;
%x0 = [0.0013;1000*2*pi]; % should be 0.0013 from w, but what for omega?
%A0 = [   ,   ;   ,   ]; %%%%%%%% TO DO
%B0 = [   ;   ]; %%%%%%%%%%%%%%%% TO DO
%[nx,nu] = size(B0);
r = [0.06; 5000*2*pi]; % reference state (maybe needs to be different)
x0 = [0.2;1000*2*pi]; % initial state (low width and high frequency does not need control)

%%% Check Controlbaility of Linear and Nonlinear Systems %%%
% If rank(ctrb(A,B)) == nx then the system is controlable
% rank(ctrb(A,B)); % need to formulate it such that C is part of Ax + Bu
% rank(ctrb(A0,B0)); % need to formulate it such that C is part of Ax + Bu

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

% Q
%Q = 1*eye(nx);
Q = 1*eye(2);

%% Terminal Constraints and Cost for Linear System

%%%Can't use DLQR because it assumes wrong cost function%%%
%[K,P,~] = dlqr(A0,B0,Q,[]);
%K_lqr = -K;

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
K = value(Y)*value(O)^-1;

%%% Constructing Terminal Set (to make Recursively Feasible)%%%
X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);
CA_set = Polyhedron(U_set.A*K_lqr,U_set.b);
IA_set = CA_set&X_set;
model = LTISystem('A',A0+B0*K_lqr);
INVset = model.invariantSet('X',IA_set); % invariant set with CA set
M_N = INVset.A;
b_N = INVset.b;

%% Linear MPC - Yalmip
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
objective = 0;
constraints = [];
for i = 1:N
    %objective = objective + x(:,i)'*Q*x(:,i)+u(:,i)'*R*u(:,i);
    objective = objective + (x(:,i)-r)'*Q*(x(:,i)-r);
    constraints = [constraints, x(:,i+1)==A0*x(:,i)+B0*u(:,i)];
    constraints = [constraints, xmin<=x(:,i)<=xmax, umin<=u(:,i)<=umax];
end
objective = objective + (x(:,N+1)-r)'*P*(x(:,N+1)-r);
constraints = [constraints, M_N*x(:,N+1)<=b_N];
MPC_Linear = optimizer(constraints,objective,[],x(:,1),u);

k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic
    Uk = MPC_Linear(xk(:,k));
    tk(:,k) = toc;
    uk(:,k) = Uk(1:nu);
    % Use discrete nonlinear model to update state
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end
data(1).x = xk;
data(1).u = uk;
data(1).t = tk;

%% Nonlinear MPC - YALMIP (solver:fmincon)
nx = 2;
nu = 1;
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
objective = 0;
constraints = [];
for i = 1:N
    %objective = objective + x(:,i)'*Q*x(:,i)+u(:,i)'*R*u(:,i);
    objective = objective + (x(:,i)-r)'*Q*(x(:,i)-r);
    % Use Nonlinear Dynamics
    A = [(1+(4/3)*(kappa*Ts*j_BS*(rho1(x(:,i),w_marg))))+(Const(1,:)/x(1,i)),  0; (-(rho2(x(:,i))*Ts)/(zeta))+(Const(1,:)/x(1,i)), (1-Ts/tau_E)+(Const(2,:)/x(2,i))];
    B = [-(kappa*Ts*eta_CD/w_dep)*rho3(x(:,i), w_dep); 0];
    
    constraints = [constraints, x(:,i+1)==A*x(:,i)+B*u(:,i)];
    constraints = [constraints, xmin<=x(:,i)<=xmax, umin<=u(:,i)<=umax];
end
%objective = objective + (x(:,N+1)-r)'*P*(x(:,N+1)-r);
%constraints = [constraints, M_N*x(:,N+1)<=b_N];

objective = objective + (x(:,N+1)-r)'*Q*(x(:,N+1)-r);
constraints = [constraints, xmin<=x(:,N+1)<=xmax];
MPC_Nonlinear = optimizer(constraints,objective,[],x(:,1),u);
% Check the solver of "MPC_Nonlinear": "Solver: FMINCON-STANDARD"

k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic
    Uk = MPC_Nonlinear(xk(:,k));
    tk(:,k) = toc;
    uk(:,k) = Uk(1:nu);
    % Use discrete nonlinear model
    A = [(1+(4/3)*(kappa*Ts*j_BS*(rho1(xk(:,k),w_marg))))+(Const(1,:)/xk(1,k)),  0; (-(rho2(xk(:,k))*Ts)/(zeta))+(Const(1,:)/xk(1,k)), (1-Ts/tau_E)+(Const(2,:)/xk(2,k))];
    B = [-(kappa*Ts*eta_CD/w_dep)*rho3(xk(:,k), w_dep); 0];
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end
data(2).x = xk;
data(2).u = uk;
data(2).t = tk;

%% Nonlinear MPC - MATLAB fmincon
U_pre = zeros(nu,N);
k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic
    [Uk, fval, exitflag] = fmincon(@(u)myfun(u,xk(:,k)),U_pre,[],[],[],[],[],[],@(u)mycon(u,xk(:,k),M_N,b_N));
    tk(:,k) = toc;
    uk(:,k) = Uk(1:nu);
    U_pre = Uk;
    % Use discrete nonlinear model
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end
data(3).x = xk;
data(3).u = uk;
data(3).t = tk;

%% Plot TEST 1
figure
subplot(2,2,1)
hold on
grid on
stairs(0:k_sim,data(2).x(1,:),'r')
%legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('x1')
subplot(2,2,2)
hold on
grid on
stairs(0:k_sim,data(2).x(2,:),'r')
%legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('x2')
subplot(2,2,3)
hold on
grid on
stairs(0:k_sim-1,data(2).u,'r')
%legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('u')
subplot(2,2,4)
hold on
grid on
plot(data(2).x(1,:),data(2).x(2,:),'r')
%legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('x1')
ylabel('x2')

figure
hold on
grid on
stairs(0:k_sim-1,data(2).t,'r')
%legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('T')

%% Plot TEST 2
figure('Position', [100 100 1000 300])
subplot(1,2,1);
plot(Ts:Ts:k_sim*Ts,data(2).u/1e6,'color',"#77AC30")
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$P_{ECCD}$ [MW]','Interpreter','latex')
title("Optimal Input")

subplot(1,2,2);
hold on
yyaxis left
plot(0:Ts:k_sim*Ts,data(2).x(1,:)*100)
plot(0:Ts:k_sim*Ts,r(1)*ones(k_sim+1)*100,'--')
ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
yyaxis right
plot(0:Ts:k_sim*Ts,data(2).x(2,:)/(2*pi))
plot(0:Ts:k_sim*Ts,r(2)*ones(k_sim+1)/(2*pi),'--')
ylabel('$\omega$ [Hz]','Interpreter','latex')
hold off
xlabel('$t$ [s]','Interpreter','latex')
legend('Output','Reference')
title("Output")

%% Plot
figure
subplot(2,2,1)
hold on
grid on
stairs(0:k_sim,data(1).x(1,:),'b')
stairs(0:k_sim,data(2).x(1,:),'r')
stairs(0:k_sim,data(3).x(1,:),'g--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('x1')
subplot(2,2,2)
hold on
grid on
stairs(0:k_sim,data(1).x(2,:),'b')
stairs(0:k_sim,data(2).x(2,:),'r')
stairs(0:k_sim,data(3).x(2,:),'g--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('x2')
subplot(2,2,3)
hold on
grid on
stairs(0:k_sim-1,data(1).u,'b')
stairs(0:k_sim-1,data(2).u,'r')
stairs(0:k_sim-1,data(3).u,'g--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('u')
subplot(2,2,4)
hold on
grid on
plot(data(1).x(1,:),data(1).x(2,:),'b')
plot(data(2).x(1,:),data(2).x(2,:),'r')
plot(data(3).x(1,:),data(3).x(2,:),'g--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('x1')
ylabel('x2')

figure
hold on
grid on
stairs(0:k_sim-1,data(1).t,'b')
stairs(0:k_sim-1,data(2).t,'r')
stairs(0:k_sim-1,data(3).t,'g--')
legend('linear MPC','Nonlinear MPC:Yalmip fmincon','Nonlinear MPC:MATLAB fmincon')
xlabel('k')
ylabel('T')