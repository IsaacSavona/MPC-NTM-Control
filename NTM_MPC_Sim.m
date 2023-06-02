%% Reset
clear all
%close all
%% Quasi-LPV MPC Model %%

%%% Model
N = 3; %prediction horizon
Ts = 0.1; %sampling time
x0 = [0;1000]; %initial state

%%% System
A = ;
B = ;

%%% Constraints
min_width = 0.06; % [m] RANDOM VALUE, 
max_width = 0.15; % [m]

%%% Constraint Polyhedrons
% state constraints (A_x * x <= b_x)
X_set = Polyhedron([-eye(nx);eye(nx)],[-xmin;xmax]);
% input constraints (A_u * u <= b_u)
U_set = Polyhedron([-eye(nu);eye(nu)],[-umin;umax]);

%%% Cost Function
Q = eye(2);

%%% Compact Formulation
Rho1 = repmat(rho1(x0(:)),1,N); % initial Rho by using current rho(k) at every predicted step
Rho2 = repmat(rho2(x0(:)),1,N);
Rho3 = repmat(rho3(x0(:)),1,N);
[Phi, Gamma] = Rho_to_PhiGamma(Rho1,Rho2,Rho3);
Omega=Q;
for j=2:(N)
  Omega=blkdiag(Omega,Q);
end
G = Gamma'*Omega*Gamma;
[W, L, c] = getWLc(A,B,xmax,xmin,umax,umin,Gamma,Phi);

%% Quasi-LPV MPC Simulation %%

%%% initialize variables
k_sim = 20;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
Uk = zeros(nu*N,k_sim);
xN = repmat(zeros(size(x0)),1,N);
Uold = ones(size(Uk));
epsilon = 1e-14;
opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym');

%%% Controller Iterations

for k = 1:k_sim % simulation loop of time samples
    for iterations = 1:10 % numerical convergence
        
            %%% Run Quadprog
            [U,~,exitflag] = quadprog(2*G,[],L,c+W*xk(:,k),[],[],[],[],[],opt);
            if exitflag ~= 1
                warning('exitflag quadprog = %d\n', exitflag)
                if exitflag == -2
                    sprintf('Optimization problem is infeasible.')
                end
            end

            %%% Store states and inputs
            Uk(:,k) = U; % stores all optimal inputs over N for each k
            uk(:,k) = U(1:nu); % pick first value of U
        
            %%% Update Rho
            xN(:,1) = xk(:,k); % Take value of the current state at k

            for i = 1:N % predict until prediction horizon N starting from state k
                xN(:,i+1) = A(Rho1(i), Rho2(i))*xN(:,i)+B(Rho3(i))*U(i); %the next predicted x is based on the optimized input sequence U
                % Update rho values with new predicted state
                Rho1(i) = rho1(xN(:,i));
                Rho2(i) = rho2(xN(:,i));
                Rho3(i) = rho3(xN(:,i));

            end
            
            [Phi, Gamma] = Rho_to_PhiGamma(Rho1,Rho2,Rho3); % Update Compact Formulation
            G = Gamma'*Omega*Gamma; % Update cost function matrix

            if (sum(abs(Uold - Uk(:,k))) < epsilon) %if the previous Uk(i-1) is less than the current Uk(i)
                disp("converged at iterations " + iterations) %display at which iteration the sufficient Rho was found
                Uold = ones(size(Uk));  %Reset Uold for next itteration
                break; %then no need to iterate further
            end
            Uold = Uk(:,k); %set what will be the previous Uk(i-1) as the value of the current Uk(i)
    end
    
    xk(:,k+1) = A(rho1, rho2)*xk(:,k)+B(rho3)*uk(:,k); % evolve state one step into future
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
legend("x1","x2","x3","x4");
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
legend("u1","u2");
% legend("First Generator Mechanical Power [W]", ...
%     "Second Generator Mechanical Power [W]", ...
%    'Interpreter','latex')

