clear;
close all;

addpath /home/sss1/Desktop/projects/dictionary/;
addpath /home/sss1/Desktop/projects/dictionary/sptoeplitz/;

num_trials = 100;
N = 1000;
s = floor(sqrt(N));
K = 1;
n = 10;
sigma = 0.1;
num_iterations = 500;
lambdas = logspace(-2, 4, 13);

load /home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/experiment3.mat;

figure;
hold all;
plot(lambdas, mean(result, 1), '*-');
% errorbar(lambdas, mean(result, 1), 2 * std(result, [], 1) ./ sqrt(num_trials));
plot(lambdas, ub, '^-');
plot(lambdas, lb, 'v-');
plot([s s], [10.^(-4) 10.^2], 'r--');
ylim([10.^(-4) 10.^2]);
notex = s/2;
notey = 10.^(mean(log10([min(lb) max(ub)])));
text(notex, notey, '$\lambda = \|R\|_{1,1}$', 'Rotation', 90, 'Interpreter', 'latex', 'FontSize', 20);
xlabel('Sparsity Parameter ($\lambda$)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Average L_2 Error', 'FontSize', 20);
set(gca, 'XScale', 'log', 'YScale', 'log');
legend({'$\widehat{X}_{\lambda}$', 'Upper bound', 'Lower bound'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northwest');
fig_root = '/home/sss1/Desktop/projects/dictionary/reconstruction_synthetic_experiments/figs/';
fig_name = [fig_root 'experiment3'];
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
saveas(gcf, [fig_name '.fig']);
saveas(gcf, [fig_name '.png']);
set(gcf, 'Units', 'Inches'); pos = get(gcf, 'Position'); set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, [fig_name '.pdf'], '-dpdf', '-r0')
% exit;
