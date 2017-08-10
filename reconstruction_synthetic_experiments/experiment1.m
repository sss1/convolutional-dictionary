clear;
close all;

addpath /home/sss1/Desktop/projects/dictionary/;
addpath /home/sss1/Desktop/projects/dictionary/sptoeplitz/;

num_trials = 100;
Ns = round(logspace(2, 4, 13)); % Sample sizes
K = 2;
n = 10;
sigma = 0.1;
num_iterations = 200;

Ss = {@(N) 5, @(N) round(N^(1/2)), @(N) round(N/10)};
Ss_names = {'$s = 5$', '$s = \lfloor N^{-1/2} \rfloor$',  '$s = \lfloor N/10 \rfloor$'}; % for plotting labels

result = zeros(num_trials, length(Ns), length(Ss)); % L_2 error in each trial
result_trivial_zero = zeros(num_trials, length(Ns), length(Ss)); % L_2 error of trivial estimator X_hat = 0
result_trivial_X = zeros(num_trials, length(Ns), length(Ss)); % L_2 error of trivial estimator X_hat = X
ub = zeros(length(Ss), length(Ns)); % upper bound
lb = zeros(length(Ss), length(Ns)); % lower bound

for N_idx = 1:length(Ns)
  N = Ns(N_idx);
  for S_idx = 1:length(Ss)
    s = Ss{S_idx}(N); % sparsity

    disp(['N: ' num2str(N) ' s: ' num2str(s)]);
    ub(S_idx, N_idx) = upper_bound(N, n, s, sigma);
    lb(S_idx, N_idx) = lower_bound(N, n, s, sigma);

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
    
      [R_hat, D_hat, ~] = learn_sparse_dictionary(X, n, K, 'lambda_D_smooth', 0, 'num_iterations', num_iterations, 'zero_edge', false);
      X_hat = multiconv(R_hat, D_hat); % estimate
    
      result(trial, N_idx, S_idx) = sum((X_hat - X_0).^2)./N;
      result_trivial_zero(trial, N_idx, S_idx) = sum(X_0.^2)./N;
      result_trivial_X(trial, N_idx, S_idx) = sum((X - X_0).^2)./N;
    end
  end
end

save('/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/experiment1.mat', 'result', 'result_trivial_zero', 'result_trivial_X', 'ub', 'lb');

figure;
set(gca,'fontsize',24);
for S_idx = 1:length(Ss)
  subplot(1, length(Ss), S_idx);
  hold all;
  errorbar(Ns, mean(result(:, :, S_idx), 1), 2 * std(result(:, :, S_idx) ./ sqrt(num_trials), [], 1));
  errorbar(Ns, mean(result_trivial_zero(:, :, S_idx), 1), 2 * std(result_trivial_zero(:, :, S_idx) ./ sqrt(num_trials), [], 1));
  errorbar(Ns, mean(result_trivial_X(:, :, S_idx), 1), 2 * std(result_trivial_X(:, :, S_idx) ./ sqrt(num_trials), [], 1));
  plot(Ns, ub(S_idx, :));
  plot(Ns, lb(S_idx, :));
  set(gca, 'XScale', 'log', 'YScale', 'log');
  title(Ss_names{S_idx});
  xlim([0.9*min(Ns) 1.1*max(Ns)]);
  ylim([0.001 1]);

  if S_idx == 1
    ylabel('Average L_2 Error');
    legend({'$\widehat{X}_s$', '$\widehat{X}_0$', '$\widehat{X}_{\infty}$', 'Upper bound', 'Lower bound'}, 'Interpreter', 'latex');
  end
  if S_idx == 2
    xlabel('Sample Size');
  end
end
fig_root = '/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/figs/';
fig_name = [fig_root 'experiment1'];
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
saveas(gcf, [fig_name '.fig']);
saveas(gcf, [fig_name '.png']);
set(gcf, 'Units', 'Inches'); pos = get(gcf, 'Position'); set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, [fig_name '.pdf'], '-dpdf', '-r0')
% exit;
