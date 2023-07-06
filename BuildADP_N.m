function OCP = BuildADP_N(sys)
% Written by D. Krishnamoorthy, Oct 2018

import casadi.*

global nx nu nd u0 lbx ubx lbu ubu

%lbu = -1e5.*ones(nu,1);
%ubu = 1e5.*ones(nu,1);

%% Build NLP solver

% empty nlp
w   = {};
w0  = [];
lbw = [];
ubw = [];

J   = 0;

g   = {}; % Inequality constraint
lbg = [];
ubg = [];

discrete = [];

x_plot  = {};
u_plot  = {};

x_init  = MX.sym('x_init',nx);
r = MX.sym('r',nx);
%d = MX.sym('d',nd);
%%

Xk = x_init;
Q = [100 0; 0 0];
for k = 0:sys.N-1
    
    Uk  = MX.sym(['U_' num2str(k)],nu);
    w   = {w{:},Uk};
    lbw = [lbw;lbu];
    ubw = [ubw;ubu];
    w0  = [w0;u0];
    u_plot   = {u_plot{:},Uk};
    %discrete = [discrete;ones(nu,1)];
    
    J = J + (Xk-r)'*Q*(Xk-r);
    %J = J + sys.J(Xk,Uk,d);
    %J = J + sys.J(Xk,Uk);

    %Xk = sys.F(Xk,Uk,d);
    Xk = sys.F(Xk,Uk);
    
    %g   = {g{:},sys.D*Uk-sys.bu}; % inequality constraint
%     g   = {g{:},Xk-ubx}; % inequality constraint
%     %lbg = [lbg;-Inf*ones(numel(sys.bu),1)];
%     lbg = [lbg;-Inf*ones(nx,1)];
%     %ubg = [ubg;zeros(numel(sys.bu),1)];
%     ubg = [ubg;zeros(nx,1)];
%     g   = {g{:},-Xk+lbx}; 
%     lbg = [lbg;-Inf*ones(nx,1)];
%     ubg = [ubg;zeros(nx,1)];
    
end

%V = (Xk-r)'*sys.P*(Xk-r);  % Terminal penalty

%%

opts = struct('warn_initial_bounds',false, ...
    'print_time',false, ...
    'ipopt',struct('print_level',1) ...
    );

w =  vertcat(w{:});  % [u0,u1,...,u_{N-1}]
G = vertcat(g{:});

%nlp  = struct('x',vertcat(w),'p',vertcat(x_init,d),'f',J+V,'g',G);
nlp  = struct('x',vertcat(w),'p',vertcat(x_init,r),'f',J,'g',G);

%if sys.discrete
    %solver = nlpsol('solver', 'bonmin', nlp, struct('discrete', discrete));
%else
    solver  = nlpsol('solver','ipopt',nlp,opts);
%end

%%

OCP.x0 = w0;
OCP.lbx = lbw;
OCP.ubx = ubw;
OCP.lbg = lbg;
OCP.ubg = ubg;
OCP.nlp = nlp;
OCP.solver = solver;

