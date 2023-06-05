%% Reset
clear all

%% Quasi-LPV MPC Model %%

%%% Model
N = 3;    % prediction horizon
Ts = 0.1; % sampling time
nx = 2;   % dimensions of state vector
nu = 1;   % dimensions of input vector
x0 = [0;1000]; % initial state (low width and high frequency does not need control)

%%% System
%A = ;
%B = ;

%%% Constraints
min_width = 0.06; % [m] RANDOM VALUE based on 2.4cm deposition width + error in dep. pos.
max_width = 0.15; % [m]
min_freq = 100;   % [Hz]
max_freq = 5000;  % [Hz]
xmin = [min_width;min_freq]; % minimum on state vector
xmax = [max_width;max_freq]; % maximum on state vector

min_power = 0;    % [W]
max_power = 2e6;  % [W]
umin = min_power; % minimum on input vector
umax = max_power; % maximum on input vector

%%% Constraint Polyhedrons
% state constraints (A_x * x <= b_x)
X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
% input constraints (A_u * u <= b_u)
U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);

%%% Cost Function
Q = eye(nx);

%%% Compact Formulation
Rho1 = repmat(rho1(x0(:)),1,N); % initial Rho by using current rho(k) at every predicted step
Rho2 = repmat(rho2(x0(:)),1,N);
Rho3 = repmat(rho3(x0(:)),1,N);
[Phi, Gamma] = Rho_to_PhiGamma(Rho1,Rho2,Rho3);
Omega = Q; % compact notation of Q is a (N*nx)x(N*nx) matrix with Q on the (block)diagonal
for j=2:(N)
  Omega = blkdiag(Omega,Q);
end
G = Gamma'*Omega*Gamma;
[W, L, c] = getWLc(xmax,xmin,umax,umin,Gamma,Phi); % Contraint matrices


%% Quasi-LPV MPC Simulation %%

%%% initialize variables
k_sim = 20; % number of simulation time steps
i_sim = 10; % max allowed number of iterations to reach numericla convergenge
xk = [x0 zeros(nx,k_sim)]; % states [w(k),ω(k)] at every time step k=0 ... k=k_sim (size=(nx) x (k_sim+1))
uk = zeros(nu,k_sim);      % input vectors [P_ECCD(k)] at every time step k=1 ... k=k_sim (size=(nu) x (k_sim))
Uk = zeros(nu*N,k_sim);    % all N predicted inputs at each k_sim time steps (size=(N) x ((nu) x (k_sim)))
xN = repmat(zeros(size(x0)),1,N); % N predicted inputs at current k only (size=(nx) x (N))
Uold = ones(size(Uk));     % Uk from previous (numerical convergence) iteration
epsilon = 1e-14;           % maximum allowed numerical error (|Uk-Uold)|)
opt =  optimoptions('quadprog','Display','off'); % create optimization options
%warning('off','optim:quadprog:HessianNotSym');   % warn if things go bad??

%%% Controller Iterations

for k = 1:k_sim              % simulation loop over time samples
    for iterations = 1:i_sim % loop until numerically convergenced (or failed)
        
            %%% Run Quadprog
            [U,~,exitflag] = quadprog(2*G,[],L,c+W*xk(:,k),[],[],[],[],[],opt); % optimize inputs U for prediction horizon given system and constraints
            if exitflag ~= 1 % if quadprog failed, give a warning
                warning('exitflag quadprog = %d\n', exitflag)
                if exitflag == -2
                    sprintf('Optimization problem is infeasible.')
                end
            end

            %%% Store states and inputs
            Uk(:,k) = U;       % stores all optimal inputs over prediction horizon for current k
            uk(:,k) = U(1:nu); % pick first value of U, the optimal input at the current time step
        
            %%% Update Rho
            xN(:,1) = xk(:,k); % take value of the current state at k (x_{0|k}=x_{k})

            for i = 1:N % predict until prediction horizon N starting from state k
                xN(:,i+1) = A(Rho1(i), Rho2(i))*xN(:,i)+B(Rho3(i))*U(i); % the next predicted x is based on the optimized input sequence U
                Rho1(i) = rho1(xN(:,i)); % update rho values with new predicted state
                Rho2(i) = rho2(xN(:,i));
                Rho3(i) = rho3(xN(:,i));
            end
            
            [Phi, Gamma] = Rho_to_PhiGamma(Rho1,Rho2,Rho3); % update Compact Formulation
            G = Gamma'*Omega*Gamma;                         % update cost function matrix

            if (sum(abs(Uold - Uk(:,k))) < epsilon)           % numerical convergence if change in Uk compared to previous iteration (Uold) is small
                disp("converged at iterations " + iterations) % display at which iteration the sufficient Rho was found
                break;                                        % if numerically converged, then no need to iterate further
            end
            Uold = Uk(:,k); % set the new previous Uk(i-1) as the current Uk(i)
    end
    
    xk(:,k+1) = A(rho1, rho2)*xk(:,k)+B(rho3)*uk(:,k); % evolve state one time step
end


%% Plot 4.1.1

%%% Plots the states and input trajectories over time%%%
figure
subplot(1,2,1)
xkPlot = xk';
hold on
stairs(0:k_sim,xk')
hold off
xlabel('$k$','Interpreter','latex')
ylabel('$x$ Angle and Frequency Deviation','Interpreter','latex')
% title('Constrained quasi-LPV MPC State Trajectory')
%axis([0 k_sim -5 10])
legend("$\mathrm{w}$ [m]","$\omega$ [Hz]");
% legend("First Generator Angle Deviation [rad]", ...
%     "First Generator Frequency Deviation [rad/s]", ...
%     "Second Generator Angle Deviation [rad]", ...
%     "Second Generator Frequency Deviation [rad/s]", ...
%    'Interpreter','latex')
subplot(1,2,2)
hold on
stairs(0:k_sim-1,uk')
hold off
xlabel('$k$','Interpreter','latex')
ylabel('$u$','Interpreter','latex')
sgtitle('Constrained quasi-LPV MPC State and Input Trajectory')
%axis([0 k_sim -1 1])
legend("$P_{ECCD}$ [W]");
% legend("First Generator Mechanical Power [W]", ...
%     "Second Generator Mechanical Power [W]", ...
%    'Interpreter','latex')

