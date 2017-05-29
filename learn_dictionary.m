% Decomposes a (column vector) signal X into the sum over j = 1:k of
% convolutions of sparse `position' vectors R_j with short, smooth `feature' vectors D_j
% Inputs:
%   X (real, N X 1): signal to decompose into R and D such that sum_j R_j * D_j
%   n (integer, scalar): desired length of features
%   K (integer, scalar): desired number of features
% Outputs:
%   R (non-negative, real, sparse, (N - n + 1) X K): matrix whose jth column
%     encodes the positions at which feature j occurs
%   D (smooth, real, n X K): matrix whose j^th column encodes the j^th feature
function [R, D, reconstruction_error] = learn_dictionary(X, n, K)

  N = length(X);
  '[N n K]'
  [N n K]

  % Initialize each D_j as a uniformly random unit vector (with 0's at each end)
  D = normc([zeros(1, K); normrnd(0, 1, [(n - 2) K]); zeros(1, K)]);

  % Initialize R to retain low frequency components of abs(X).
  % We then add random noise (to break symmetry between columns, and then
  % enforce non-negativity.
  smoothing_kernel = normpdf((1:n) - mean(1:n), 0, n^(1/2));
  R = repmat(conv(abs(X), smoothing_kernel, 'valid'), [1 K]);
  R = max(normrnd(R, 1), 0);
  '[length(smoothing_kernel) size(R, 1) length(X)]'
  [length(smoothing_kernel) size(R, 1) length(X)]

  lambda_D_smooth = 5; % weight of smoothness penalty on D
  lambda_R_sparse = 1; % weight of sparsity penalty on R
  num_iterations = 5000; % total number of gradient descent iterations
  gamma_0 = 0.1; % initial gradient descent step size

  reconstruction_error = zeros(num_iterations, 1);
  for iteration = 1:num_iterations
    iteration
    gamma = gamma_0 .* (iteration.^(-1/2));


    D = D - gamma * grad_D(D, R, X, lambda_D_smooth);
    D = normc(D); % ensure that columns of D are unit vectors
    R = R - gamma * grad_R(D, R, X, lambda_R_sparse);
    R = max(R, 0); % ensure that R is non-negative

    reconstruction_error(iteration) = norm(sum(multiconv(R, D), 2) - X);

  end

end

% Helper function: calculates the gradient of the objective with respect to D
function D_grad = grad_D(D, R, X, lambda_D_smooth)

  [n, K] = size(D);

  % Gradient of penalty for roughness in D
  smoothness_grad =[zeros(1, K); ...
                    (2.*D(2:(n - 1), :)) - D(1:(n - 2), :) - D(3:n, :); ...
                    zeros(1, K)];

  % Gradient of error in reconstructing X
  cost_grad = zeros(size(D));
  % 'size(sum(multiconv(R, D), 2))'
  % size(sum(multiconv(R, D), 2))
  % 'size(X)'
  % size(X)
  w = 2 .* (sum(multiconv(R, D), 2) - X)';
  for k = 1:K
    cost_grad(:, k) = w * convmtx(R(:, k), n);
  end

  D_grad = cost_grad + smoothness_grad;

  D_grad(1, :) = 0; D_grad(n, :) = 0; % fix end points of D to zero

end

% Helper function: calculates the gradient of the objective with respect to R
function R_grad = grad_R(D, R, X, lambda_R_sparse)

  [n, K] = size(D);

  % gradient for L1 norm of R
  R_grad = zeros(size(R));
  sparsity_grad = lambda_R_sparse;

  % Gradient of error in reconstructing X
  cost_grad = zeros(size(R));
  w = 2 .* (sum(multiconv(R, D), 2) - X)';
  for k = 1:K
    cost_grad(:, k) = w * convmtx(D(:, k), size(R, 1));
  end

  R_grad = cost_grad + sparsity_grad;

end
