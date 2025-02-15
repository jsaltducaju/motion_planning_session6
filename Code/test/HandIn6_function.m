function [out] = Handin_function(state_i, state_f)
% Compute Dubins curves using optimization with the tool CasADi
% using direct collocation for discretization of the continuous-time
% motion equations. - Function version for testing

%Check input size
a = size(state_i);
b = size(state_f);
if a(1) ~= 3 || a(2) ~= 1
   error('HandIn6:InputSize','state_i must be a matrix of size 3x1')
end
if b(1) ~= 3 || b(2) ~= 1
   error('HandIn6:InputSize','state_f must be a matrix of size 3x1')
end


%Check input variable type
if ~isa(state_i(1),'numeric') || ~isa(state_i(2),'numeric') || ~isa(state_i(3),'numeric')
    error('HandIn6:InputVariableType','state_i elements must be numbers.')
end
if ~isa(state_f(1),'numeric') || ~isa(state_f(2),'numeric') || ~isa(state_f(3),'numeric')
    error('HandIn6:InputVariableType','state_f elements must be numbers.')
end

%Input constraint
u_max = 1; %u_min = -u_max

% Parameters for collocation
N = 75; % Number of elements
nx = 3; % Degree of state vector
Nc = 3; % Degree of interpolation polynomials

% Formulate the optimization problem for minimum path length using CasADi
addpath('../casadi-osx-matlabR2015a-v3.5.5')
import casadi.*

% Use the opti interface in CasADi
opti = casadi.Opti();

% Define optimization variables and motion equations
x = MX.sym('x',nx);
u = MX.sym('u');
f = Function('f',{x, u}, {cos(x(3)), sin(x(3)), u});
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
ds = T_opt; %ds = T*v (v=1)
 
%Matlab Dubins 
dubConnObj = dubinsConnection;
[pathSegObj, pathCosts] = connect(dubConnObj,state_i',state_f');


%Check solution
out = abs(pathCosts - T_opt);

end