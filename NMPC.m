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
N = 2;                % prediction horizon
Ts = 0.1;             % controller sampling time
%r = [0.00125; 5000*2*pi]; % reference state (maybe needs to be different)
%r1 = [0.0835; 5000*2*pi]; % reference state (maybe needs to be different)
%x0 = [0.06;1000*2*pi]; % initial state (low width and high frequency does not need control)
x0 = [0.04; 1000*2*pi];
%r1 = [0.04;5000*2*pi];
%x0 = [0.10; 1000*2*pi];
r1 = [0.086;5000*2*pi];
r2 = [0.10;5000*2*pi];
r3 = [0.03;5000*2*pi];
%x0 = r;
%Q = [300 0; 0 1e-8];        % no frequency error to reference taken into account!
Q = [300 0; 0 0];        
%Q = 100;
R = 0.0001;
%R = 0;
R1 = 0.001;

nx = 2; % dimensions of state
nu = 1; % dimensions of input


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
ustar = sdpvar(nu);
K = sdpvar(nu+nx,1);
objective = 0;
constraints = [];
%load eq;
C=[1 0];
% SS reference 1
A01 = A(rho1(r1,w_marg), rho2(r1), kappa,Ts,j_BS,zeta,tau_E);
B01 = B(rho3(r1,w_dep), kappa,Ts,eta_CD,w_dep);
K1=inv([A01-eye(2) B01; C 0])*[Const; r1(1)];
%K1(3)=0.516;
% SS reference 2
A02 = A(rho1(r2,w_marg), rho2(r2), kappa,Ts,j_BS,zeta,tau_E);
B02 = B(rho3(r2,w_dep), kappa,Ts,eta_CD,w_dep);
K2=inv([A02-eye(2) B02; C 0])*[Const; r2(1)];
% SS referenece 3
A03 = A(rho1(r3,w_marg), rho2(r3), kappa,Ts,j_BS,zeta,tau_E);
B03 = B(rho3(r3,w_dep), kappa,Ts,eta_CD,w_dep);
K3=inv([A03-eye(2) B03; C 0])*[Const; r3(1)];
r = sdpvar(nx,1);


for i = 1:N
    %objective = objective + (x(:,i)-K(1:2))'*Q*(x(:,i)-K(1:2)) + (u(:,i)-K(3))'*R*(u(:,i)-K(3)) + (u(:,i)-ustar)'*R1*(u(:,i)-ustar);
    objective = objective + (x(:,i)-r)'*Q*(x(:,i)-r);
    %objective = objective + (C*x(:,i)-K(1))'*Q*(C*x(:,i)-K(1)) + (u(:,i)-K(3))'*R*(u(:,i)-K(3));% + (u(:,i)-u(:,i-1))'*R1*(u(:,i)-u(:,i-1));
    % Use Nonlinear Dynamics
    Amat = A(rho1(x(:,i),w_marg), rho2(x(:,i)), kappa,Ts,j_BS,zeta,tau_E);
    Bmat = B(rho3(x(:,i),w_dep), kappa,Ts,eta_CD,w_dep);
    constraints = [constraints,x(:,i+1)==Amat*x(:,i)+Bmat*u(:,i)+Const, ... % the state equation should always be met
                   umin<=u(:,i)<=umax, x(:,i)<=xmax, xmin(2)<=x(2,i)'];%, ... % input and state constraints
                   %diff(u(:,i))]; % constraint of change in input to prevent 
    
end

% Apply terminal constraints
%objective = objective + (x(:,N+1)-K(1:2))'*Q*(x(:,N+1)-K(1:2));
objective = objective + (x(:,N+1)-r)'*Q*(x(:,N+1)-r);
constraints = [constraints, xmin(2)<=x(2,N+1)<=xmax(2), x(1,N+1)<=xmax(1)];
% Define optimizer
options = sdpsettings('usex0',0,'solver','fmincon','verbose',0);
optimoptions('fmincon','Algorithm','sqp');
%MPC_Nonlinear = optimizer(constraints,objective,options,{x(:,1),ustar,K},u);
MPC_Nonlinear = optimizer(constraints,objective,options,{x(:,1),ustar,r},u);
% Check the solver of "MPC_Nonlinear": "Solver: FMINCON-STANDARD"


%% Simulation
k_sim = 250; % number of simulation time steps
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
ukstar = uk(:,1);
%U0 = sdpvar(nu,N);
%assign(u,2);
%Kk = K1;
rk = r1;
for k = 1:k_sim
    k
    %disp([r1(1),Ref(k)])
    if k == 125
        %Kk = K2;
        rk = r2;
        %disp([r2(1),Ref(k)])
    end
    if k == 175
        %Kk = K3;
        rk = r3;
        %disp([r3(1),Ref(k)])
    end
    
    % Every sampling time (starting at 1), adjust controller input
    tic
    %[Uk,~,~,~,~] = MPC_Nonlinear({xk(:,k),ukstar,Kk});
    [Uk,~,~,~,~] = MPC_Nonlinear({xk(:,k),ukstar,rk});
    tk(k) = toc;
    %assign(u,[Uk(:,2:N),Kk(3)])
    %assign(U0,[Uk(1,2:N);K(3)])
    uk(:,k) = Uk(1:nu);
    ukstar = uk(:,k);
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

subplot(2,2,1);
plot(Ts:Ts:k_sim*Ts,uk(:),'color',"#77AC30")
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$P_{ECCD}$ [MW]','Interpreter','latex')
title(sprintf("Optimal input with $N =$ %d, $T_s =$ %0.1e s",N,Ts),'Interpreter','latex')

subplot(2,2,3);
hold on
yyaxis left
plot(0:Ts:k_sim*Ts,xk(1,:)*100)
%plot(0:k_sim,xk(1,:)*100)
Ref = [r1(1)*ones(1,124)*100 r2(1)*ones(1,49)*100 r3(1)*ones(1,78)*100];
%plot(0:k_sim,Ref,'--')
stairs(0:Ts:k_sim*Ts,Ref,'--')
% plot(0:Ts:124*Ts,r1(1)*ones(125)*100,'--')
% xline(124*Ts);
% plot(124*Ts:Ts:174*Ts,r2(1)*ones(50+1)*100,'--')
% xline(174*Ts);
% plot(174*Ts:Ts:k_sim*Ts,r3(1)*ones(75+2)*100,'--')
ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
yyaxis right
plot(0:Ts:k_sim*Ts,xk(2,:)/(2*pi))
ylabel('$\omega$ [Hz]','Interpreter','latex')
hold off
xlabel('$t$ [s]','Interpreter','latex')
legend('Output','Reference')
title(sprintf("Output with $\\mathrm{w}_0 =$ %0.3f cm, $\\omega_0 =$ %0.0f Hz",x0(1)*100,x0(2)/2/pi),'Interpreter','latex')

subplot(2,2,[2,4]);
plot(1:k_sim,data(1).t,'color',"#7E2F8E")
xlabel('$k$ [samples]','Interpreter','latex')
ylabel('$T_{comp}$ [seconds]','Interpreter','latex')
title("Time to Compute Optimal Solution with SQP",'Interpreter','latex')
