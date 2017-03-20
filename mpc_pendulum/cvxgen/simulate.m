function simulate()

    load('datafile.mat');

    fprintf('\nSolving MPC pendulum problem loop for N = %d and solver cvxgen', N)
    
    whos

    x0 = vec(x0);
    params.A = full(A);
    params.B = full(B);
    params.Q = full(Q);
    params.QN = full(QN);
    params.R = full(R);
    params.x_0 = vec(full(x0));
    params.x_max = vec(full(xmax));
    params.u_max = vec(full(umax));
    [m,n] = size(B);

    % Initialize time vector
    time = zeros(nsim,1);

    % Initialize number of iterations vector
    niter = zeros(nsim,1);

    if     N == 10
        cd T10;
    elseif N == 20
        cd T20;
    end

    nprob = 10;

    % Setup settings
    settings.verbose = 0;

    %make_csolve;

    for i = 1:nsim

        % Solve
        tic;
        vars = csolve(params, settings);
        time(i) = toc;

        % Apply first control input to the plant
        x0 = A*x0 + B*vars.u_0;

    end

    cd ..
    save('cvxgen_results.mat', 'time');

end
