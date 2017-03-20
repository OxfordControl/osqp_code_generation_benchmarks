clear

function [ time_vec ] = cvxgen_pendulum_T10()
% Run mpc_pendulum example with T = 10 using CVXGEN

if ~exist(strcat('csolve.', mexext), 'file')
    % Generate code with
    cvxgen(333447049216)
end

% Load data
data = load('.mat');

% Setup parameters struct
params.A = full(data.F);
params.B = full(data.D);
params.Q = eye(size(data.F, 2));
params.R = data.mu';
params.QN = data.mu';
params.x0 = data.mu';
params.xmax = data.mu';
params.umax = data.mu';

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




    if 0
      m = 3;
      n = 6;
      T = 20;
      dt = .1;
      Q = randn(n,n); Q = Q'*Q;
      R = randn(m,m); R = R'*R;
      A = eye + dt*randn(n,n);
      B = dt*randn(n,m);
      umax = 20;
      x0 = randn(n,1)/dt;
    else
      load mpc_data
    end

    params.A = A;
    params.B = B;
    params.Q = Q;
    params.R = R;
    params.u_max = umax;
    params.x_0 = x0;

    tic
    [vars, status] = csolve(params);
    toc

end
