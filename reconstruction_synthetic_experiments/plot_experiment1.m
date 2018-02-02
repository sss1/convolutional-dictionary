clear;
close all;

num_trials = 100;
Ns = round(logspace(2, 4, 13)); % Sample sizes

Ss = {@(N) 5, @(N) round(N^(1/2)), @(N) round(N/10)};
Ss_names = {'$\|R\|_{1,1} = 5$', '$\|R\|_{1,1} = \lfloor N^{-1/2} \rfloor$', '$\|R\|_{1,1} = \lfloor N/10 \rfloor$'}; % for plotting labels

load /home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/experiment1.mat;

figure;
for S_idx = 1:length(Ss)
  subplot(1, length(Ss), S_idx);
  hold all;
  plot(Ns, mean(result(:, :, S_idx), 1), '*-');
  plot(Ns, mean(result_trivial_zero(:, :, S_idx), 1), 'o-');
  plot(Ns, mean(result_trivial_X(:, :, S_idx), 1), 'x-');
  % errorbar(Ns, mean(result(:, :, S_idx), 1), 2 * std(result(:, :, S_idx) ./ sqrt(num_trials), [], 1));
  % errorbar(Ns, mean(result_trivial_zero(:, :, S_idx), 1), 2 * std(result_trivial_zero(:, :, S_idx) ./ sqrt(num_trials), [], 1));
  % errorbar(Ns, mean(result_trivial_X(:, :, S_idx), 1), 2 * std(result_trivial_X(:, :, S_idx) ./ sqrt(num_trials), [], 1));
  plot(Ns, ub(S_idx, :), '^-');
  plot(Ns, lb(S_idx, :), 'v-');
  set(gca, 'XScale', 'log', 'YScale', 'log');
  title(Ss_names{S_idx}, 'Interpreter', 'latex', 'FontSize', 16);
  xlim([0.9*min(Ns) 1.1*max(Ns)]);
  ylim([0.00001 1]);

  if S_idx == 1
    ylabel('Average L_2 Error', 'FontSize', 16);
  elseif S_idx == 2
    xlabel('Sample Size', 'FontSize', 16);
  elseif S_idx == 3
    legend({'$\widehat{X}_{\lambda}$', '$\widehat{X}_0$', '$\widehat{X}_{\infty}$', 'Upper bound', 'Lower bound'}, ...
           'Location', 'southeast', ...
           'Interpreter', 'latex', ...
           'FontSize', 16);
  end
end
fig_root = '/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/figs/';
fig_name = [fig_root 'experiment1'];
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
saveas(gcf, [fig_name '.fig']);
saveas(gcf, [fig_name '.png']);
set(gcf, 'Units', 'Inches'); pos = get(gcf, 'Position'); set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, [fig_name '.pdf'], '-dpdf', '-r0')
% exit;
