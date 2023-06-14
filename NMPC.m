%% Reset
clear all


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
Cw = 1e3;      % [#] UNKNOWN!!
tau_A0 = 3e-6; % [s] Alfvèn time
tau_w = 0.188; % [s] resistive wall time
omega0 = 2*pi*420; % [rad/s] equilibrium frequency

kappa = 16*mu0*Lq*rs^2/(0.82*tau_r*B_pol*pi);
zeta = m*Cw*tau_A0^2*tau_w*a^3;


%% Controller and simulation parameters
N = 20;                % prediction horizon
Ts = 0.01;             % controller sampling time
r = [0.00125; 5000*2*pi]; % reference state (maybe needs to be different)
x0 = [0.10;1000*2*pi]; % initial state (low width and high frequency does not need control)
Q = [1 0; 0 0];        % no frequency error to reference taken into account!

nx = 2; % dimensions of state
nu = 1; % dimensions of input

k_sim = 100; % number of simulation time steps


%% Constraints
min_width = 0.00; % [m] RANDOM VALUE based on 2.4cm deposition width + error in dep. pos.
max_width = 0.15; % [m]
min_freq = 100*2*pi;  % [rad/s]
max_freq = 5000*2*pi; % [rad/s]
xmin = [min_width;min_freq]; % minimum on state vector
xmax = [max_width;max_freq]; % maximum on state vector

min_power = 0;    % [MW]
max_power = 2;    % [MW]
umin = min_power; % minimum on input vector
umax = max_power; % maximum on input vector


%% Nonlinear MPC - YALMIP (solver:fmincon)
Const = [-4/3*(kappa*Ts*j_BS*w_sat)/(w_sat^2+w_marg^2); omega0*Ts/tau_E0];
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
objective = 0;
constraints = [];
for i = 1:N
    objective = objective + (x(:,i)-r)'*Q*(x(:,i)-r);
    % Use Nonlinear Dynamics
    Amat = A(rho1(x(:,i),w_marg), rho2(x(:,i)), kappa,Ts,j_BS,zeta,tau_E);
    Bmat = B(rho3(x(:,i),w_dep), kappa,Ts,eta_CD,w_dep);
    constraints = [constraints,x(:,i+1)==Amat*x(:,i)+Bmat*u(:,i)+Const*Ts, ... % the state equation should always be met
                   umin<=u(:,i)<=umax, x(:,i)<=xmax, xmin(2)<=x(2,i)'];%, ... % input and state constraints
                   %diff(u(:,i))]; % constraint of change in input to prevent 
    
end

% Apply terminal constraints
objective = objective + (x(:,N+1)-r)'*Q*(x(:,N+1)-r);
constraints = [constraints, xmin(2)<=x(2,N+1)<=xmax(2), x(1,N+1)<=xmax(1)];
% Define optimizer
options = sdpsettings('verbose',0);
MPC_Nonlinear = optimizer(constraints,objective,options,x(:,1),u);
% Check the solver of "MPC_Nonlinear": "Solver: FMINCON-STANDARD"


%% Simulation
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    
    % Every sampling time (starting at 1), adjust controller input
    tic
    [Uk,~,~,~,~] = MPC_Nonlinear(xk(:,k));
    uk(:,k) = Uk(1:nu);
    disp(uk(k))
    tk(k) = toc;

    % Use discrete nonlinear model
    xk(:,k+1) = A(rho1(xk(:,k),w_marg), rho2(xk(:,k)), kappa,Ts,j_BS,zeta,tau_E)*xk(:,k) ...
                + B(rho3(xk(:,k),w_dep), kappa,Ts,eta_CD,w_dep)*uk(:,k) ...
                + Const;

    if xk(1,k+1) <=0
        disp('island has been destroyed')
        xk(1,k+1) = 0;
        xk(:,k+2:k_sim+1) = nan;
        uk(:,k+2:k_sim) = nan;
        tk(k+2:k_sim) = nan;
        break
    end

end
data(1).x = xk;
data(1).u = uk;
data(1).t = tk;


%% Plot output and optimal input
figure('Position', [100 100 1000 300])

subplot(1,2,1);
plot(Ts:Ts:k_sim*Ts,uk(:),'color',"#77AC30")
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$P_{ECCD}$ [MW]','Interpreter','latex')
title(sprintf("Optimal input with $N =$ %d, $T_s =$ %0.1e s",N,Ts),'Interpreter','latex')

subplot(1,2,2);
hold on
yyaxis left
plot(0:Ts:k_sim*Ts,xk(1,:)*100)
plot(0:Ts:k_sim*Ts,r(1)*ones(k_sim+1)*100,'--')
ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
yyaxis right
plot(0:Ts:k_sim*Ts,xk(2,:)/(2*pi))
ylabel('$\omega$ [Hz]','Interpreter','latex')
hold off
xlabel('$t$ [s]','Interpreter','latex')
legend('Output','Reference')
title(sprintf("Output with $\\mathrm{w}_0 =$ %0.3f cm, $\\omega_0 =$ %0.0f Hz",x0(1)*100,x0(2)/2/pi),'Interpreter','latex')