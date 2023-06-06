%% Comments
%%% Latex file corrections
% w_marg =/= w_crit, is the value correct than? --> w_marg is in fact the
% value below which the neoclassical island drive is limited (??)

%%% Sanity checks
% w_sat is indeed a stable equilibrium for island width
% w_crit = 1.25mm =/= (w_marg = 2cm) !!
% w_crit scales with w_marg^2/w_sat
% w_crit is unstable equilibrium
% ws/(ws^2+wm^2) =  w/(w^2+wm^2) --> ws*w^2 + ws*wm^2 = w*(ws^2+wm^2) -->
% w^2 - w*(ws^2+wm^2)/ws + wm^2 = 0 --> 
% (w-(ws^2+wm^2)/(2*ws))^2 = ((ws^2+wm^2)/(2*ws))^2 - wm^2 -->
% wc = (ws^2+wm^2)/(2*ws) + sqrt(((ws^2+wm^2)/(2*ws))^2 - wm^2)
% w>wc: island is growing and slowing down
% w<wc: island is vanishing and slowing down
% w~ws: w is behaving as expected and omega is blowing up quickly
% If omega<0 --> ODE solution unstable and explodes after a short wiggle
% Without constants the critical island size is 0 as there is no 
% (passively) stabilizing (constant) term. Saturation behavior, strangly
% enough, is also gone completely without the constant term.


%% Reset
clear all;

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

Delta_0 = -16*mu0*Lq/(B_pol*pi) * 4/3 * w_sat/(w_sat^2+w_marg^2)*j_BS; % NC drive

kappa = 16*mu0*Lq*rs^2/(0.82*tau_r*B_pol*pi);
zeta = m*Cw*tau_A0^2*tau_w*a^3;


%% Simulation parameters
Ts = 1e-3;          % [s] sampling time
k_sim = 3000;       % [#] number of simulation time steps
wi = 0.10;          % [m] initial island size
omegai = 2*pi*1000; % [rad/s] initial frequency
w = [wi,zeros(1,k_sim)];         % island width array as function of time
omega = [omegai,zeros(1,k_sim)]; % frequency array as function of time
P_ec = zeros(1,k_sim);           % ECCD power input
U_ECCD = 0e6;       % [W] constant control input

% %% Simulation: Rutherford & La Haye equations
% for k = 1:k_sim
%     Delta_BS = (16*mu0*Lq)/(B_pol*pi) * 4/3 * w(k)/(w(k)^2+w_marg^2)*j_BS;
%     Delta_CD = -(16*mu0*Lq)/(B_pol*pi)*(eta_CD*P_ec(k))/(w_dep)*(0.25+0.24*w(k)/w_dep)/(1+1.5*w(k)/w_dep+0.43*(w(k)/w_dep)^2+0.64*(w(k)/w_dep)^3);
%     w(k+1) = w(k) + Ts*rs/(0.82*tau_r) * (rs*Delta_0 + rs*Delta_BS + rs*Delta_CD);
%     omega(k+1) = omega(k) + Ts*(omega0/tau_E0 - omega(k)/tau_E - 1/(m*Cw*tau_A0^2*tau_w*omega(k))*(w(k)/a)^3);
%     if w(k+1) <= 0 || omega(k+1) <= 0
%         w(k+1) = 0;
%         omega(k+1) = 0;
%         break
%     end
% end

% %% Simulation: Simplified Rutherford & La Haye equations
% for k = 1:k_sim
%     rho1 = 1/(w(k)^2+w_marg^2);
%     rho2 = (0.25+0.24*w(k)/w_dep)/(1+1.5*w(k)/w_dep+0.43*(w(k)/w_dep)^2+0.64*(w(k)/w_dep)^3);
%     rho3 = w(k)^2/omega(k);
% 
%     w(k+1) = w(k)*(1 + 4/3*kappa*Ts*j_BS*rho1) ...
%            - P_ec(k)*(eta_CD*kappa*Ts)/(w_dep)*rho2 ...
%            - 4/3*(kappa*Ts*j_BS*w_sat)/(w_sat^2+w_marg^2);
%     omega(k+1) = omega(k)*(1-Ts/tau_E) - w(k)*Ts/zeta*rho3 + Ts*omega0/tau_E0;
%     if w(k+1) <= 0 || omega(k+1) <= 0
%         w(k+1) = 0;
%         omega(k+1) = 0;
%         break
%     end
% end

%% Simulation: State space notation
xc = [[wi;omegai],zeros(2,k_sim)];
u = ones(1,k_sim)*U_ECCD;
C = [- 4/3*(kappa*Ts*j_BS*w_sat)/(w_sat^2+w_marg^2) ; ...
     Ts*omega0/tau_E0];

for k = 1:k_sim
    
    rho1 = 1/(xc(1,k)^2+w_marg^2);
    rho2 = xc(1,k)^2/xc(2,k);
    rho3 = (0.25+0.24*xc(1,k)/w_dep)/(1+1.5*xc(1,k)/w_dep+0.43*(xc(1,k)/w_dep)^2+0.64*(xc(1,k)/w_dep)^3);

    A = [(1 + 4/3*kappa*Ts*j_BS*rho1)  0 ; ...
         -Ts/zeta*rho2  (1-Ts/tau_E)];
    B = [-(eta_CD*kappa*Ts)/(w_dep)*rho3 ; ...
         0];

    xc(:,k+1) = A*xc(:,k) + B*u(k) + C;

    if xc(1,k+1) <= 0
        xc(1,k+1) = 0;
        xc(:,k+2:k_sim+1) = NaN;
        break
    end
    if xc(2,k+1) <= 0
        xc(2,k+1) = 0;
        xc(:,k+2:k_sim+1) = NaN;
        break
    end
end

%% Simulation: State space without constants
x = [[wi;omegai],zeros(2,k_sim)];
u = ones(1,k_sim)*U_ECCD;

for k = 1:k_sim
    
    rho1 = 1/(x(1,k)^2+w_marg^2);
    rho2 = x(1,k)^2/x(2,k);
    rho3 = (0.25+0.24*x(1,k)/w_dep)/(1+1.5*x(1,k)/w_dep+0.43*(x(1,k)/w_dep)^2+0.64*(x(1,k)/w_dep)^3);

    A = [(1 + 4/3*kappa*Ts*j_BS*rho1)  0 ; ...
         -Ts/zeta*rho2  (1-Ts/tau_E)];
    B = [-(eta_CD*kappa*Ts)/(w_dep)*rho3 ; ...
         0];

    x(:,k+1) = A*x(:,k) + B*u(k);

    if x(1,k+1) <= 0
        x(1,k+1) = 0;
        x(:,k+2:k_sim+1) = NaN;
        break
    end
    if x(2,k+1) <= 0
        x(2,k+1) = 0;
        x(:,k+2:k_sim+1) = NaN;
        break
    end
end

% %% Visualisation w/omega
% figure
% subplot(1,1,1)
% hold on
% yyaxis left
% plot(0:Ts:k_sim*Ts,w*100)
% ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
% yyaxis right
% plot(0:Ts:k_sim*Ts,omega/(2*pi))
% ylabel('$\omega$ [Hz]','Interpreter','latex')
% hold off
% xlabel('$t$ [s]','Interpreter','latex')

% %% Visualisation x
% figure
% subplot(1,1,1)
% hold on
% yyaxis left
% plot(0:Ts:k_sim*Ts,x(1,:)*100)
% ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
% yyaxis right
% plot(0:Ts:k_sim*Ts,x(2,:)/(2*pi))
% ylabel('$\omega$ [Hz]','Interpreter','latex')
% hold off
% xlabel('$t$ [s]','Interpreter','latex')

%% Visualisation x/xc
figure
subplot(1,1,1)
hold on
yyaxis left
plot(0:Ts:k_sim*Ts,x(1,:)*100)
plot(0:Ts:k_sim*Ts,xc(1,:)*100)
ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
yyaxis right
plot(0:Ts:k_sim*Ts,x(2,:)/(2*pi))
plot(0:Ts:k_sim*Ts,xc(2,:)/(2*pi))
ylabel('$\omega$ [Hz]','Interpreter','latex')
hold off
xlabel('$t$ [s]','Interpreter','latex')
legend('Without constants','With constants')
