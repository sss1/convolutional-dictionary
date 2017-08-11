clear;
close all;

addpath /home/sss1/Desktop/projects/dictionary/;
addpath /home/sss1/Desktop/projects/dictionary/sptoeplitz/;

num_trials = 100;
N = 1000;
s = round(sqrt(N));
K = 1;
n = 10;
sigma = 0.1;
num_iterations = 500;
lambdas = logspace(0.5, 3.5, 13);

result = zeros(num_trials, length(lambdas)); % L_2 error in each trial
ub = zeros(1, length(lambdas)); % upper bound
lb = zeros(1, length(lambdas)); % lower bound

for lambda_idx = 1:length(lambdas)
  lambda = lambdas(lambda_idx);

  disp(['lambda: ' num2str(lambda)]);
  if lambda >= s % upper bound only holds assuming lambda >= s
    ub(lambda_idx) = upper_bound(N, n, lambda, sigma);
  end
  lb(lambda_idx) = lower_bound(N, n, s, sigma);

  parfor trial = 1:num_trials

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
  
    [R_hat, D_hat, ~] = learn_constrained_sparse_dictionary(X, n, K, 'lambda_D_smooth', 0, 'lambda_R_sparse', lambda, 'num_iterations', num_iterations, 'zero_edge', false);
    X_hat = multiconv(R_hat, D_hat); % estimate
  
    result(trial, lambda_idx) = sum((X_hat - X_0).^2)./N;
  end
end

save('/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/experiment3.mat', 'result', 'ub', 'lb');
plot_experiment3;
