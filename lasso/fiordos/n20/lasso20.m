%function [ time_vec ] = lasso20()
%LASSO20 Run portfolio example with n = 20, m = 200
% with FiOrdOs


% Load data
data = load('datafilen20.mat');
[m, n] = size(data.Ad);
In = eye(n);


%% Solve with OSQP

gamma = data.gammas(1);

P = blkdiag(zeros(n), speye(m), zeros(n));
q = [zeros(n+m, 1); gamma*ones(n, 1)];
A = [data.Ad, -eye(m), zeros(m, n);
     -In, zeros(n, m), -In;
     In, zeros(n, m), -In];
l = [data.bd'; zeros(2*n, 1)];
u = [data.bd'; Inf*ones(2*n, 1)];


% Initialize OSQP object
osqpSolver = osqp;
osqpSolver.setup(P, q, A, l, u);

% Solve
tic
res_osqp = osqpSolver.solve();
toc



%% Function for updating linear cost
q_func = @(gamma, n, m) [zeros(n+m, 1); gamma*ones(n, 1)];

% Simple set
X = SimpleSet(EssRn(2*n + m));
op = OptProb('H', diag([zeros(n,1); ones(m, 1); zeros(n, 1)]), 'g', 'param', 'X', X, ...
             'Ae', full([data.Ad, -eye(m), zeros(m, n)]), 'be', data.bd', ...
             'Ai', [-In, zeros(n, m), -In; In, zeros(n, m), -In], 'bi', zeros(2*n, 1));
s = Solver(op, 'approach', 'primal-dual', 'algo', 'fgm');

% Dual approach settings
s.setSettings('approach', 'inlineH', 1);
s.setSettings('approach', 'apprMaxit', 50000);

% Algorithm settings
s.setSettings('algo', 'stopg', true, 'stopgEps', 1e-1);

% Generate code
fiordos_dir = 'code';
s.generateCode('prefix', 'demo_', 'outDir', fiordos_dir, 'forceOverwrite', true);

cd(fiordos_dir)
demo_mex_make();

% Init parameters and settings
mparams = struct();
msetgs = struct();
res = struct();
res.x = zeros(2*n + m, 1);
res.la = zeros(2*n + m, 1);

% Number of problems to be solved
n_prob = length(data.gammas);

% Initialize time vector
time_vec = zeros(n_prob, 1);


% Run loop
for i = 1:n_prob
    gamma = data.gammas(i);

    tic
    mparams.g = q_func(gamma, n, m);        % Update linear cost
    msetgs.approach.apprInitX = res.x;      % Warm starting
    msetgs.approach.apprInitLa = res.la;
    
    % Solve
    res = demo_mex(mparams, msetgs);
    
    time_vec(i) = toc;

    if res.exitflag ~= 2
        fprintf('problem not solved\n');
    end

end
cd ..

%end
