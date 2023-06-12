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

Ts = 1e-2;
k_sim = 1000;

Const = [-4/3*(kappa*Ts*j_BS*w_sat)/(w_sat^2+w_marg^2); Ts*omega0/tau_E0];

%%
x = [[0.00125;1000*2*pi] zeros(2,k_sim)];

for i = 1:k_sim
    A = [(1+(4/3)*(kappa*Ts*j_BS*(rho1(x(:,i),w_marg))))+(Const(1,:)/x(1,i)),  0; (-(rho2(x(:,i))*Ts)/(zeta)), (1-Ts/tau_E)+(Const(2,:)/x(2,i))];
    %B = [-(kappa*Ts*eta_CD/w_dep)*rho3(x(:,i), w_dep); 0];
    
    x(:,i+1) = A*x(:,i);
end

%% Plot w to dwdt
plot((x(1,2:k_sim+1)+x(1,1:k_sim))/2,(x(1,2:k_sim+1)-x(1,1:k_sim))/Ts)


%% Plot

hold on
yyaxis left
plot(0:Ts:k_sim*Ts,x(1,:)*100)
%plot(0:Ts:k_sim*Ts,r(1)*ones(k_sim+1)*100,'--')
ylabel('$\mathrm{w}$ [cm]','Interpreter','latex')
yyaxis right
plot(0:Ts:k_sim*Ts,x(2,:)/(2*pi))
%plot(0:Ts:k_sim*Ts,r(2)*ones(k_sim+1)/(2*pi),'--')
ylabel('$\omega$ [Hz]','Interpreter','latex')
hold off
xlabel('$t$ [s]','Interpreter','latex')
legend('Output','Reference')
title("Output")

%% Linearizing
syms w omega
R = ((1+(4/3)*(kappa*Ts*j_BS*(1/(w^2+w_marg^2)))) - ...
    (kappa*Ts*eta_CD/w_dep)*((0.25 + 0.24*w/w_dep)/(1 + 1.5*w/w_dep + 0.43*(w/w_dep)^2 + 0.64*(w/w_dep)^3)))*w + ...
    Const(1,:);
H = ((-((w^2/omega)*Ts)/(zeta))*w + (1-Ts/tau_E))*omega + Const(2,:);
% A = [(1+(4/3)*(kappa*Ts*j_BS*(rho1(x(:,i),w_marg))))+(Const(1,:)/x(1,i)),  0; (-(rho2(x(:,i))*Ts)/(zeta)), (1-Ts/tau_E)+(Const(2,:)/x(2,i))];
% B = [-(kappa*Ts*eta_CD/w_dep)*rho3(x(:,i), w_dep); 0];
diff(R,w)
diff(R,omega)
diff(H,w)
diff(H,omega)