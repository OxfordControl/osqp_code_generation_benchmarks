% Generate and solve the problem with OSQP

clear
clc

q_func = @(gamma, mu, k) [-mu'/(2*gamma); zeros(k, 1)];

% Number of assets
n = 200;

% Load data
filename = ['../datafilen',num2str(n),'.mat'];
load(filename);

gamma = gammas(10);
k = size(F, 2); 

% OSQP data
P = blkdiag(D, speye(k));
q = q_func(gamma, mu, k);
A = [ones(1, n), zeros(1, k); F', -speye(k); speye(n), zeros(n, k)];
l = [1; zeros(k, 1); zeros(n, 1)];
u = [1; zeros(k, 1); ones(n, 1)];

% Initialize OSQP object
osqpSolver = osqp;
osqpSolver.setup(P, q, A, l, u, 'eps_rel', 1e-3, 'eps_abs', 1e-3, ...
                 'rho', 1e-1, 'sigma', 1e-5);

% Solve
tic
res_osqp = osqpSolver.solve();
toc


%% Solve the problem with Fiordos (dense)

% % Simple set
% X1 = SimpleSet(EssRnplus(n));
% op1 = OptProb('H', full(D + F*F'), 'g', 'param', 'X', X1, 'Ae', ones(1, n), 'be', 1.0);
% s1 = Solver(op1, 'approach', 'dual', 'algoOuter', 'fgm');
% 
% % Algorithm settings
% s1.setSettings('algoOuter', 'stopg', true, 'stopgEps', 1e-3);
% s1.setSettings('algoOuter', 'init', 0.0, 'maxit', 1000);
% 
% % Approach settings
% s1.setSettings('approach', 'inlineA', 1);
% 
% % Generate code
% fiordos_dir = 'fiordos_dense';
% s1.generateCode('prefix', 'demo_', 'outDir', fiordos_dir, 'forceOverwrite', true);
% 
% % Set parameters
% mparams = struct();
% mparams.g = mu' / (2*gamma);
% 
% cd(fiordos_dir)
% demo_mex_make();
% 
% tic
% res1_fiordos = demo_mex(mparams);
% toc
% cd ..

%% Solve the problem with Fiordos (sparse)

% Simple set
X2 = SimpleSet(EssBox(n, 'l', zeros(n,1), 'u', ones(n,1)), EssRn(k));
op2 = OptProb('H', full(P), 'g', 'param', 'X', X2, ...
              'Ae', [ones(1, n), zeros(1, k); full(F'), -eye(k)], 'be', [1.0; zeros(k,1)]);
s2 = Solver(op2, 'approach', 'dual', 'algoOuter', 'fgm');

% Algorithm settings
s2.setSettings('algoOuter', 'stopg', true, 'stopgEps', 1e-3);
s2.setSettings('algoOuter', 'init', zeros(k+1, 1), 'maxit', 10000);

% Approach settings
s2.setSettings('approach', 'inlineA', 1);

% Generate code
fiordos_dir = 'fiordos_sparse';
s2.generateCode('prefix', 'demo_', 'outDir', fiordos_dir, 'forceOverwrite', true);

% Set parameters
mparams = struct();
mparams.g = q_func(gamma, mu, k);

cd(fiordos_dir)
demo_mex_make();

tic
res2_fiordos = demo_mex(mparams);
toc
cd ..
