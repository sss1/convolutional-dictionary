clear;
close all;

addpath /home/sss1/Desktop/projects/dictionary/;
addpath /home/sss1/Desktop/projects/dictionary/sptoeplitz/;

num_trials = 300;
Ns = round(logspace(2, 4, 13)); % Sample sizes
K = 1;
n = 30;
sigma = 1;
p = 2; % Number of finite moments
num_iterations = 500;

Ss = {@(N) 5, @(N) round(N^(1/2)), @(N) round(N/10)};
Ss_names = {'$s = 5$', '$s = \lfloor N^{-1/2} \rfloor$',  '$s = \lfloor N/10 \rfloor$'}; % for plotting labels

result = zeros(num_trials, length(Ns), length(Ss)); % L_2 error in each trial
result_trivial_zero = zeros(num_trials, length(Ns), length(Ss)); % L_2 error of trivial estimator X_hat = 0
result_trivial_X = zeros(num_trials, length(Ns), length(Ss)); % L_2 error of trivial estimator X_hat = X
ub = zeros(length(Ss), length(Ns)); % upper bound
% lb = zeros(length(Ss), length(Ns)); % lower bound

for N_idx = 1:length(Ns)
  N = Ns(N_idx);
  for S_idx = 1:length(Ss)
    s = Ss{S_idx}(N); % sparsity

    disp(['N: ' num2str(N) ' s: ' num2str(s)]);
    ub(S_idx, N_idx) = moment_upper_bound(N, n, s, sigma, p);
    % lb(S_idx, N_idx) = lower_bound(N, n, s, sigma);

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
      signs = ((rand(N,1)<.5)*2 - 1);
      epsilon = signs .* gprnd(1/p, sigma, sigma*p, N, 1); % Pareto noise (with 2-\eps moments) with random signs
      X = X_0 + epsilon; % data signal
    
      [R_hat, D_hat, ~] = general_dictionary(X, n, K, 'lambda_R_sparse', s, 'num_iterations', num_iterations);
      X_hat = multiconv(R_hat, D_hat); % estimate
    
      result(trial, N_idx, S_idx) = sum((X_hat - X_0).^2)./N;
      result_trivial_zero(trial, N_idx, S_idx) = sum(X_0.^2)./N;
      result_trivial_X(trial, N_idx, S_idx) = sum((X - X_0).^2)./N;
    end
  end
end

save('/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/experiment4.mat', 'result', 'result_trivial_zero', 'result_trivial_X', 'ub');%, 'lb');
plot_experiment4;
