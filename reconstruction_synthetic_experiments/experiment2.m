clear;
close all;

addpath /home/sss1/Desktop/projects/dictionary/;
addpath /home/sss1/Desktop/projects/dictionary/sptoeplitz/;

num_trials = 100;
N = 1000;
s = 10;
K = 1;
ns = round(logspace(0.5, 2, 10));
sigma = 0.1;
num_iterations = 500;

result = zeros(num_trials, length(ns)); % L_2 error in each trial
result_trivial_zero = zeros(num_trials, length(ns)); % L_2 error of trivial estimator X_hat = 0
result_trivial_X = zeros(num_trials, length(ns)); % L_2 error of trivial estimator X_hat = X
ub = zeros(1, length(ns)); % upper bound
lb = zeros(1, length(ns)); % lower bound

for n_idx = 1:length(ns)
  n = ns(n_idx);

  disp(['n: ' num2str(n)]);
  ub(n_idx) = upper_bound(N, n, s, sigma);
  lb(n_idx) = lower_bound(N, n, s, sigma);

  parfor trial = 1:num_trials

    % D = normc(ones(n, K));
    D = normc(normrnd(0, 1, n, K));
    R = zeros(N - n + 1, K);
  
    for i = 1:s
      for k = 1:K
        idx = randi(N - n + 1);
        R(idx, k) = R(idx, k) + 1;
      end
    end
  
    X_0 = multiconv(R, D); % estimand
    X = normrnd(X_0, sigma); % data signal
  
    [R_hat, D_hat, ~] = learn_constrained_sparse_dictionary(X, n, K, 'lambda_D_smooth', 0, 'lambda_R_sparse', s, 'num_iterations', num_iterations, 'zero_edge', false);
    X_hat = multiconv(R_hat, D_hat); % estimate
  
    result(trial, n_idx) = sum((X_hat - X_0).^2)./N;
    result_trivial_zero(trial, n_idx) = sum(X_0.^2)./N;
    result_trivial_X(trial, n_idx) = sum((X - X_0).^2)./N;
  end
end

save('/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/experiment2.mat', 'result', 'result_trivial_zero', 'result_trivial_X', 'ub', 'lb');
plot_experiment2;
