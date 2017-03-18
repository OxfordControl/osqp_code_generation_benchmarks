function filename = pendulum_data( N )
%CREATE_DATA Summary of this function goes here
%   N -- prediction horizon

filename = ['data/pendulum_N', num2str(N), '.mat'];

if exist(filename, 'file') == 2
    disp(['File ', filename, ' has already been generated.']);
    return;
end

% Model data
A = [1., 0.02, 0., 0.;
	0., 1., 0., 0.;
    0., 0., 1.0049, 0.02;
    0., 0., 0.4913, 1.0049];
B = [0.0002, 0.02, -0.0005, -0.0501]';
[nx, nu] = size(B);

% MPC data
Q = diag([1., 0.3, 0.3, 0.1]);
R = 0.1;
QN = dare(A,B,Q,R);
% QN = Q;

% Initial state
x_init = [0., 0., 0.15, 0.]';
x_lim = [0.5, 1.0, 0.2, 0.5]';

% Input constraints
D = 1.0;
dl = -5.0;
du = 5.0;

% State constraints
E = eye(nx);
el = -x_lim;
eu = x_lim;

% Terminal state constraints
F = E;
fl = el;
fu = eu;


%% MPC matrices

% Objective
Hx = kron(eye(N), Q);
Hu = kron(eye(N), R);
H = blkdiag(Hx, QN, Hu);

% Dynamics
Mx = kron(eye(N+1), -eye(nx)) + kron(diag(ones(N,1), -1), A);
Mu = kron([zeros(1,N); eye(N)], B);
M = sparse([Mx, Mu]);

% Constraints
Gx = kron(eye(N), E);
Gu = kron(eye(N), D);
G = sparse(blkdiag(Gx, F, Gu));
gl = [repmat(el, N, 1); fl; repmat(dl, N, 1)];
gu = [repmat(eu, N, 1); fu; repmat(du, N, 1)];


%% Define the QP data

b = @(x, n, N) [-x; zeros(N*n, 1)];

P = sparse(H);
q = zeros((N+1)*nx + N*nu, 1);
A = sparse([M; G]);
l = [b(x_init, nx, N); gl];
u = [b(x_init, nx, N); gu];

% Save the data
save(filename, 'P', 'q', 'A', 'M', 'l', 'u', 'gl', 'gu', 'b', 'nx', 'nu', 'N', 'x_init');

return

end

