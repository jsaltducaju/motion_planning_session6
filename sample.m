addpath('../casadi-osx-matlabR2015a-v3.5.5')

% Compute the motion primitives using optimization with the tool CasADi
% using direct collocation for discretization of the continuous-time
% motion equations.
% Parameters for collocation
N = 75; % Number of elements
nx = 3; % Degree of state vector
Nc = 3; % Degree of interpolation polynomials
x_vec = 2;
y_vec = 2;
th_vec = -pi;
v = 1;
L = 1;
u_max = 1;
u_min = -1;
state_i = [0;0;0];
% Formulate the optimization problem for minimum path length using CasADi
import casadi.*
for i = 1:length(x_vec)

 % Use the opti interface in CasADi
 opti = casadi.Opti();
 state_f = [x_vec(i) y_vec(i) th_vec(i)]';
 % Define optimization variables and motion equations
 x = MX.sym('x',nx);
 u = MX.sym('u');
 f = Function('f',{x, u}, {cos(x(3)), sin(x(3)), tan(u)});
 X = opti.variable(nx,N+1);
 pos_x = X(1,:);
 pos_y = X(2,:);
 ang_th = X(3,:);
 U = opti.variable(N,1);
 T = opti.variable(1);
 % Set the element length (with final time T unknown, and thus an
 % optimization variable)
 dt = T/N;
 % Set initial guess values of variables
 opti.set_initial(T,0.1);
 opti.set_initial(U,0.0*ones(N,1));

 % Define collocation parameters
 tau = collocation_points(Nc,'radau');
 [C,~] = collocation_interpolators(tau);
 % Formulate collocation constraints
 for k = 1:N % Loop over elements
 Xc = opti.variable(nx,Nc);
 X_kc = [X(:,k) Xc];
 for j = 1:Nc
 % Make sure that the motion equations are satisfied at
 % all collocation points
 [f_1, f_2, f_3] = f(Xc(:,j),U(k));
 opti.subject_to(X_kc*C{j+1}' == dt*[f_1; f_2; f_3]);
 end
 % Continuity constraints for states between elements
 opti.subject_to(X_kc(:,Nc+1) == X(:,k+1));
 end
 % Input constraints
 for k = 1:N
 opti.subject_to(-u_max <= U(k) <= u_max);
 end
 % Initial and terminal constraints
 opti.subject_to(T >= 0.001);
 opti.subject_to(X(:,1) == state_i);
 opti.subject_to(X(:,end) == state_f);
 % Formulate the cost function
 alpha = 1e-2;
 opti.minimize(T);
 % Choose solver ipopt and solve the problem
 opti.solver('ipopt',struct('expand',true),struct('tol',1e-8));
 sol = opti.solve();
 % Extract solution trajectories and store them in mprim variable
 pos_x_opt = sol.value(pos_x);
 pos_y_opt = sol.value(pos_y);
 ang_th_opt = sol.value(ang_th);
 u_opt = sol.value(U);
 T_opt = sol.value(T);
 mprim{i}.x = pos_x_opt;
 mprim{i}.y = pos_y_opt;
 mprim{i}.th = ang_th_opt;
 mprim{i}.u = u_opt;
 mprim{i}.T = T_opt;
 mprim{i}.ds = T_opt*v;
 end