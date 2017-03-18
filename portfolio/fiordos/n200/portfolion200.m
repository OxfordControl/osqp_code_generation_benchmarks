function [ time_vec ] = portfolion200()
%PORTFOLION100CVXGEN Run portfolio example with n = 200, k = 20
% with CVXGEN


% Load data
data = load('datafilen200.mat');
[n, k] = size(data.F);

% Function for updating linear cost
q_func = @(gamma, mu, k) [-mu'/(2*gamma); zeros(k, 1)];

% Simple set
X = SimpleSet(EssBox(n, 'l', zeros(n, 1), 'u', ones(n,1)), EssRn(k));
op = OptProb('H', blkdiag(full(data.D), eye(k)), 'g', 'param', 'X', X, ...
             'Ae', [ones(1, n), zeros(1, k); full(data.F'), -eye(k)], 'be', [1.0; zeros(k,1)]);
s = Solver(op, 'approach', 'dual', 'algoOuter', 'fgm');

% Dual approach settings
s.setSettings('approach', 'inlineA', 1);

% Algorithm settings
s.setSettings('algoOuter', 'stopg', true, 'stopgEps', 1e-4);
s.setSettings('algoOuter', 'init', zeros(k+1, 1), 'maxit', 50000);

% Generate code
fiordos_dir = 'code';
s.generateCode('prefix', 'demo_', 'outDir', fiordos_dir, 'forceOverwrite', true);

cd(fiordos_dir)
demo_mex_make();

% Init parameters and settings
mparams = struct();
msetgs = struct();
res = struct();
res.la = zeros(k+1, 1);

% Number of problems to be solved
n_prob = length(data.gammas);

% Initialize time vector
time_vec = zeros(n_prob, 1);


% Run loop
for i = 1:n_prob
    gamma = data.gammas(i);

    tic
    mparams.g = q_func(gamma, data.mu, k);   % Update linear cost
    msetgs.algoOuter.init = res.la;     % Warm starting
    res = demo_mex(mparams, msetgs);
    time_vec(i) = toc;

    if res.exitflag ~= 2
        fprintf('problem not solved\n');
    end

end
cd ..




end
