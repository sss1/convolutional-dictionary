%clear;
close all;

load /home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/experiment2.mat;

figure;
hold all;
plot(ns, mean(result, 1), '*-');
plot(ns, mean(result_trivial_zero, 1), 'o-');
plot(ns, mean(result_trivial_X, 1), 'x-');
% errorbar(ns, mean(result, 1), 2 * std(result, [], 1) ./ sqrt(num_trials));
% errorbar(ns, mean(result_trivial_zero, 1), 2 * std(result_trivial_zero, [], 1) ./ sqrt(num_trials));
% errorbar(ns, mean(result_trivial_X, 1), 2 * std(result_trivial_X, [], 1) ./ sqrt(num_trials));
plot(ns, ub, '^-');
plot(ns, lb, 'v-');
xlabel('Dictionary Element Length ($n$)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Average L_2 Error', 'FontSize', 20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([0.9*min(ns) 1.1*max(ns)]);
ylim([0.9*min(lb) 0.1]);
legend({'$\widehat{X}_s$', '$\widehat{X}_0$', '$\widehat{X}_{\infty}$', 'Upper bound', 'Lower bound'}, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 20);
fig_root = '/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/figs/';
fig_name = [fig_root 'experiment2'];
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
saveas(gcf, [fig_name '.fig']);
saveas(gcf, [fig_name '.png']);
set(gcf, 'Units', 'Inches'); pos = get(gcf, 'Position'); set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, [fig_name '.pdf'], '-dpdf', '-r0')
% exit;
