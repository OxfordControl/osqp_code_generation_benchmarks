function [ time_vec ] = portfolion100cvxgen()
%PORTFOLION100CVXGEN Run portfolio example with n = 100, m = 10
% with CVXGEN


# Generate code with
% cvxgen(526211878912)

% Load data
data = load('datafilen100.mat');

% Setup parameters struct
params.F = full(data.F);
params.D = full(data.D);
params.I = eye(size(data.F, 2));
params.mu = data.mu';

% Setup settings
settings.verbose = 0;

% Number of problems to be solved
n_prob = length(data.gammas);

% Initialize time vector
time_vec = zeros(n_prob, 1);


% Run loop
for i = 1:n_prob
    params.gamma = data.gammas(i);

    tic
    [vars, status] = csolve(params, settings);
    time_vec(i) = toc;

    if status.converged ~= 1
        fprintf('problem not solved\n');
    end

%     [vars_cvx, status_cvx] = cvxsolve(params, settings);

%     pause(0.1);


end





end
