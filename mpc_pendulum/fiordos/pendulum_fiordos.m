% Generate and solve the problem with OSQP

clear
clc

% Prediction horizon
N = 3;

% Generate and load data
if ~(exist('data', 'dir') == 7)
   mkdir data 
end

filename = pendulum_data(N);
load(filename);

% Initialize OSQP object
osqpSolver = osqp;
osqpSolver.setup(P, q, A, l, u, 'eps_rel', 1e-3, 'eps_abs', 1e-3, ...
                 'rho', 1e-1, 'sigma', 1e-3, 'alpha', 1.95);

% Solve
tic
res_osqp = osqpSolver.solve();
toc


%% Solve the problem with Fiordos

[m, n] = size(M);

% Simple set
X = SimpleSet(EssBox(n, 'l', gl, 'u', gu));
op = OptProb('H', full(P), 'g', q, 'X', X, 'Ae', M, 'be', 'param');
s = Solver(op, 'approach', 'dual', 'algoOuter', 'fgm');

% Algorithm settings
s.setSettings('algoOuter', 'stopg', true, 'stopgEps', 1e-3);
s.setSettings('algoOuter', 'init', zeros(m, 1), 'maxit', 10000);

% Approach settings
s.setSettings('approach', 'inlineA', 1);

% Generate code
fiordos_dir = 'fiordos_code';
s.generateCode('prefix', 'demo_', 'outDir', fiordos_dir, 'forceOverwrite', true);

% Set parameters
mparams = struct();
mparams.be = b(x_init, nx, N);

cd(fiordos_dir)
demo_mex_make();

tic
res_fiordos = demo_mex(mparams);
toc
cd ..
