%% Reset
clear all

%% Physics parameters
j_BS = 73e3;   % [A/m^2] bootstrap current density
w_dep = 0.024; % [m] deposition width
w_marg = 0.02; % [m] critical width
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


%% Quasi-LPV MPC Model %%

%%% Model
N = 20;    % prediction horizon
Ts = 0.01; % sampling time
nx = 2;   % dimensions of state vector
nu = 1;   % dimensions of input vector
x0 = [0.2;1000*2*pi]; % initial state (low width and high frequency does not need control)
%x0 = [2;1000*2*pi]; % initial state (low width and high frequency does not need control)

%%% System
C = [-4/3*(kappa*Ts*j_BS*w_sat)/(w_sat^2+w_marg^2); Ts*omega0/tau_E0];

%%% Constraints
% min_width = 0.06; % [m] RANDOM VALUE based on 2.4cm deposition width + error in dep. pos.
% max_width = 0.15; % [m]
% min_freq = 100*2*pi;  % [rad/s]
% max_freq = 5000*2*pi; % [rad/s]
% xmin = [min_width;min_freq]; % minimum on state vector
% xmax = [max_width;max_freq]; % maximum on state vector
% 
% min_power = 0;    % [W]
% max_power = 2e6;  % [W]
% umin = min_power; % minimum on input vector
% umax = max_power; % maximum on input vector

%%% Test Constraints
min_width = -100; % [m] RANDOM VALUE based on 2.4cm deposition width + error in dep. pos.
max_width = 1000; % [m]
min_freq = 0*2*pi;  % [rad/s]
max_freq = 50000*2*pi; % [rad/s]
xmin = [min_width;min_freq]; % minimum on state vector
xmax = [max_width;max_freq]; % maximum on state vector

min_power = 0;    % [W]
max_power = 2e6;  % [W]
umin = min_power; % minimum on input vector
umax = max_power; % maximum on input vector

% %%% Constraint Polyhedrons
% % state constraints (A_x * x <= b_x)
% X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
% % input constraints (A_u * u <= b_u)
% U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);

%%% Cost Function
%Q = 2*eye(nx);       % weights on width and freq deviation from reference
%Q = zeros(nx)
Q = [1 0; 0 0];
r = [0.06; 5000*2*pi]; % reference state
%r = [min_width; 5000*2*pi]; % reference state

%%% Compact Formulation
Rho1 = repmat(rho1(x0(:),w_marg),1,N); % compact notation of initial Rho by using current rho(k) at every predicted step
Rho2 = repmat(rho2(x0(:)),1,N);
Rho3 = repmat(rho3(x0(:),w_dep),1,N);
[Phi, Gamma, Lambda] = Rho_to_PhiGammaLambda(Rho1,Rho2,Rho3, @A,@B,C, kappa,Ts,j_BS,zeta,tau_E,eta_CD,w_dep);
Omega = Q;       % compact notation of Q is a (N*nx)x(N*nx) matrix with Q on the (block)diagonal
for j=2:N
  Omega = blkdiag(Omega,Q);
end
R = repmat(r,N,1)'; % compact notation of R
G = 2*Gamma'*Omega*Gamma;                % quadratic part of cost function (U^T G U)
F = 2*Gamma'*Omega*(Phi*x0+Lambda-R'); % linear part of cost function (F^T U)
[W, L, c] = getWLc(xmax,xmin,umax,umin,Gamma,Phi,Lambda); % contraint matrices


%% Quasi-LPV MPC Simulation %%

%%% initialize variables
k_sim = 25; % number of simulation time steps
i_sim = 10; % max allowed number of iterations to reach numerical convergence
xk = [x0 zeros(nx,k_sim)]; % states [w(k),ω(k)] at every time step k=0 ... k=k_sim (size=(nx) x (k_sim+1))
uk = zeros(nu,k_sim);      % input vectors [P_ECCD(k)] at every time step k=1 ... k=k_sim (size=(nu) x (k_sim))
Uk = zeros(nu*N,k_sim);    % all N predicted inputs at each k_sim time steps (size=(N) x ((nu) x (k_sim)))
xN = repmat(zeros(size(x0)),1,N); % N predicted inputs at current k only (size=(nx) x (N))
Uold = ones(size(Uk));     % Uk from previous (numerical convergence) iteration
epsilon = 1e-4;           % maximum allowed numerical error (|Uk-Uold)|)
opt =  optimoptions('quadprog','Display','off','MaxIterations',400); % create optimization options
%warning('off','optim:quadprog:HessianNotSym');   % warn if things go bad??

%%% Controller Iterations

for k = 1:k_sim              % simulation loop over time samples
    for iterations = 1:i_sim % loop until numerically convergenced (or failed)
       
            %%% Run Quadprog
            %[U,~,exitflag] = quadprog(G,F,L,c+W*xk(:,k),[],[],[],[],[],opt); % optimize inputs U for prediction horizon given system and constraints
            [U,~,exitflag] = quadprog(G,F,[],[],[],[],[],[],[],opt); % optimize inputs U for prediction horizon given system and constraints
            if exitflag ~= 1 % if quadprog failed, give a warning
                if exitflag == 0
                    disp('Solver stopped prematurely.')
                elseif exitflag == -2
                    sprintf('Optimization problem is infeasible.')
                else
                    warning('exitflag quadprog = %d\n', exitflag)
                end
            end

            %%% Store states and inputs
            Uk(:,k) = U;       % stores all optimal inputs over prediction horizon for current k
            uk(:,k) = U(1:nu); % pick first value of U, the optimal input at the current time step
       
            %%% Update Rho
            xN(:,1) = xk(:,k); % take value of the current state at k (x_{0|k}=x_{k})

            for i = 1:N % predict until prediction horizon N starting from state k
                xN(:,i+1) = A(Rho1(i), Rho2(i), kappa,Ts,j_BS,zeta,tau_E)*xN(:,i)+B(Rho3(i), kappa,Ts,eta_CD,w_dep)*U(i)+C; % the next predicted x is based on the optimized input sequence U
                Rho1(i) = rho1(xN(:,i),w_marg); % update rho values with new predicted state
                Rho2(i) = rho2(xN(:,i));
                Rho3(i) = rho3(xN(:,i),w_dep);
            end
           
            [Phi, Gamma, Lambda] = Rho_to_PhiGammaLambda(Rho1,Rho2,Rho3, @A,@B,C, kappa,Ts,j_BS,zeta,tau_E,eta_CD,w_dep); % update Compact Formulation
            G = 2*Gamma'*Omega*Gamma;                % update quadratic part of cost function (U^T G U)
            F = 2*Gamma'*Omega*(Phi*xk(:,k)+Lambda-R'); % update linear part of cost function (F^T U)
            [W, L, c] = getWLc(xmax,xmin,umax,umin,Gamma,Phi,Lambda); % update contraint matrices

            if (sum(abs(Uold - Uk(:,k))) < epsilon)           % numerical convergence if change in Uk compared to previous iteration (Uold) is small
                disp("converged at iterations " + iterations) % display at which iteration the sufficient Rho was found
                break;                                        % if numerically converged, then no need to iterate further
            end
            Uold = Uk(:,k); % set the new previous Uk(i-1) as the current Uk(i)
    end
   
    xk(:,k+1) = A(rho1(xk(:,k),w_marg), rho2(xk(:,k)), kappa,Ts,j_BS,zeta,tau_E)*xk(:,k)+B(rho3(xk(:,k),w_dep), kappa,Ts,eta_CD,w_dep)*uk(:,k)+C; % evolve state one time step
end


%% Plot output and input
figure('Position', [100 100 1000 300])
subplot(1,2,1);
plot(Ts:Ts:k_sim*Ts,uk(:)/1e6,'color',"#77AC30")
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$P_{ECCD}$ [MW]','Interpreter','latex')
title("Optimal Input")

subplot(1,2,2);
hold on
yyaxis left
plot(0:Ts:k_sim*Ts,xk(1,:)*100)
plot(0:Ts:k_sim*Ts,r(1)*ones(k_sim+1)*100,'--')
ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
yyaxis right
plot(0:Ts:k_sim*Ts,xk(2,:)/(2*pi))
plot(0:Ts:k_sim*Ts,r(2)*ones(k_sim+1)/(2*pi),'--')
ylabel('$\omega$ [Hz]','Interpreter','latex')
hold off
xlabel('$t$ [s]','Interpreter','latex')
legend('Output','Reference')
title("Output")