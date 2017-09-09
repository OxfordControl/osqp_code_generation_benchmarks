function [ time_vec ] = portfolion100()
%PORTFOLION100CVXGEN Run portfolio example with n = 100, m = 10
% with CVXGEN

if ~exist(strcat('csolve.', mexext), 'file')
    % Generate code with
    cvxgen(307372040192)
end

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
    [~, status] = csolve(params, settings);
    time_vec(i) = toc;

    if status.converged ~= 1
        fprintf('problem not solved\n');
    end

%     [vars_cvx, status_cvx] = cvxsolve(params, settings);

%     pause(0.1);


end





end
