clear;
close all;

sigma = 0.1; % noise variance
lambda_R_sparse = 0.1; % regularization parameter making R sparse
lambda_D_smooth = 0.5; % regularization parameter making D smooth
gamma_0_R = 0.1; % initial gradient descent step size for R
gamma_0_D = 1; % initial gradient descent step size for D

% D is a smooth feature vector
D = normc([0 2 3  2 0;...
           0 1 0 -1 0]');
[n, K] = size(D);

% R is a length R_length sparse vector taking the value 1 at R_1
% distinct, uniformly random coordinates
R_length = 100;
R_1 = 20;
R = zeros(R_length, K);
R(datasample(1:numel(R), R_1, 'Replace', false)) = 1;

% X is the Gaussian data, centered at the convolution of R and D
X = normrnd(sum(multiconv(R, D), 2), sigma);

[R_hat, D_hat, reconstruction_error] = learn_dictionary(X, n, K);
D_hat_switched = D(:, [2 1]);

% % switch the order if necessary, just for plotting
% if norm(D_hat_switched(:) - D(:)) < norm(D_hat(:) - D(:))
%   D_hat = D_hat_switched;
%   R_hat = R(:, [2 1]);
% end


X_hat = sum(multiconv(R_hat, D_hat), 2); % Reconstruction of X

figure;
subplot(5, 2, 1); plot(D(:, 1)); title('True D1');
subplot(5, 2, 2); plot(D_hat(:, 1)); title('Estimated D1');
subplot(5, 2, 3); plot(D(:, 2)); title('True D2');
subplot(5, 2, 4); plot(D_hat(:, 2)); title('Estimated D2');
subplot(5, 2, 5); stem(R(:, 1)); title('True R1');
subplot(5, 2, 6); stem(R_hat(:, 1)); title('Estimated R1');
subplot(5, 2, 7); stem(R(:, 2)); title('True R2');
subplot(5, 2, 8); stem(R_hat(:, 2)); title('Estimated R2');

subplot(5, 2, 9); hold all; title('True X');
plot(X); plot(sum(multiconv(R, D), 2)); plot(X_hat);
legend('Noisy', 'Pure', 'Estimated');
subplot(5, 2, 10); plot(X_hat); title('Estimated X');
