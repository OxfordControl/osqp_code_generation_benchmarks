% Run all CVXGEN tests


% Number of assets
% n_vec = [20, 30, 50, 80, 100, 120];
n_vec = [50, 80, 100, 120];

n_probs = length(n_vec);

% Get current folder
cur_dir = pwd;

% Initialize time vector
avg_vec = zeros(n_probs, 1);
std_vec = zeros(n_probs, 1);
median_vec = zeros(n_probs, 1);
max_vec = zeros(n_probs, 1);


% Simulate system for all problems
for i = 1:n_probs
    cd(sprintf('n%i', n_vec(i)))

    % Construct function handle
    fh = str2func(sprintf('portfolion%i', n_vec(i)));

    % Call simulation function
    time_temp = fh();

    % Compute statistics
    avg_vec(i) = mean(time_temp);
    std_vec(i) = std(time_temp);
    median_vec(i) = median(time_temp);
    max_vec(i) = max(time_temp);

    % Go back to original directory
    cd(cur_dir)


end


% Store vectors to file
save('cvxgen_results.mat', 'avg_vec', 'std_vec', 'median_vec', 'max_vec');
