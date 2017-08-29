function [R, D, reconstruction_error] ...
              = learn_constrained_sparse_dictionary(X, n, K, varargin)
%
% Decomposes a (column vector) signal X into the sum over j = 1:k of
% convolutions of sparse `position' vectors R_j with short, smooth `feature' vectors D_j
% 
% Required Inputs:
%   X (real, N X 1):        signal to decompose into R and D (such that sum_j R_j * D_j approximates X)
%   n (integer scalar > 0): desired length of features
%   K (integer scalar > 0): desired number of features
% 
% Optional Inputs:
%   lambda_D_smooth (real scalar >= 0,    Default: 0.01) weight of smoothness penalty on D
%   lambda_R_sparse (real scalar > 0,     Default: 0.01) weight of sparsity penalty on R
%   gamma_0         (real scalar >= 0,    Default: 0.01) initial gradient descent step size
%   num_iterations  (integer scalar >= 0, Default: 100)  total number of gradient descent iterations
% 
% Outputs:
%   R (non-negative, real, sparse, (N - n + 1) X K): matrix whose jth column encodes positions at which feature j occurs
%   D (smooth, real, n X K): matrix whose j^th column encodes the j^th feature
%   opt (non-negative, real, num_iterations X 1): history of objective function

% TODO: Figure out if bottleneck computations (mostly, building the (sparse) convolution matrix) can be sped up on GPU

  % default values for optional arguments
  lambda_D_smooth = 0.01;
  lambda_R_sparse = 0.01;
  gamma_0 = 0.01;
  num_iterations = 200;
  zero_edge = true;

  % parse optional arguments
  for i=1:2:length(varargin)
     switch varargin{i}
     case 'lambda_D_smooth'
       lambda_D_smooth = varargin{i + 1};
       if lambda_D_smooth < 0
         error('Optional parameter lambda_D_smooth must be non-negative');
       end
     case 'lambda_R_sparse'
       lambda_R_sparse = varargin{i + 1};
       if lambda_R_sparse < 0
         error('Optional parameter lambda_R_sparse must be non-negative');
       end
     case 'gamma_0'
       gamma_0 = varargin{i + 1};
       if gamma_0 < 0
         error('Optional parameter gamma_0 must be non-negative');
       end
     case 'num_iterations'
       num_iterations = varargin{i + 1};
     case 'zero_edge'
       zero_edge = varargin{i + 1};
     otherwise
       error(sprintf('%s is not a valid argument name', varargin{i}));
     end
  end

  N = length(X);

  % Initialize each D_j as a uniformly random unit vector (with 0's at each end)
  D = normc([zeros(1, K); normrnd(0, 1, [(n - 2) K]); zeros(1, K)]);

  % Initialize each R_j as a non-negative IID gaussian vector
  R = max(normrnd(0, var(X)./K, [(N - n + 1), K]), 0);

  if nargout >= 3
    reconstruction_error = zeros(num_iterations, 1);
  end
  for iteration = 1:num_iterations
    gamma = gamma_0 .* (iteration.^(-1/2));

    D = D - gamma * grad_D(D, R, X, lambda_D_smooth, zero_edge);
    D = normc(D); % ensure that columns of D are unit vectors
    R = R - gamma * grad_R(D, R, X);
    R = max(R, 0); % ensure that R is non-negative
    for k = 1:K
      R(:, k) = ProjectOntoL1Ball(R(:, k), lambda_R_sparse); % ensure that ||R||_1 <= lambda_R_sparse
    end

    if nargout >= 3
      reconstruction_error(iteration) = norm(multiconv(R, D) - X) + norm(R(:), 1);
    end

  end

end

% Helper function: calculates the gradient of the objective with respect to D
function D_grad = grad_D(D, R, X, lambda_D_smooth, zero_edge)

  [n, K] = size(D);

  % Gradient of penalty for roughness in D
  smoothness_grad = [zeros(1, K); ...
                     (2.*D(2:(n - 1), :)) - D(1:(n - 2), :) - D(3:n, :); ...
                     zeros(1, K)];

  % Gradient of error in reconstructing X
  cost_grad = zeros(size(D));
  w = 2 .* (multiconv(R, D) - X)';
  for k = 1:K
    % For a column vector c, sptoeplitz([c; zeros(n - 1,1)], zeros(n,1))
    % is a sparse equivalent of convmtx(c, n).
    cost_grad(:, k) = w * sptoeplitz([R(:, k); zeros(n - 1,1)], [R(1, k); zeros(n - 1,1)]);
  end

  D_grad = cost_grad + lambda_D_smooth .* smoothness_grad;
  if zero_edge
    D_grad(1, :) = 0; D_grad(n, :) = 0; % fix end points of D to zero
  end

end

% Helper function: calculates the gradient of the objective with respect to R
function R_grad = grad_R(D, R, X)

  [n, K] = size(D);

  % gradient for L1 norm of R
  R_grad = zeros(size(R));

  % Gradient of error in reconstructing X
  cost_grad = zeros(size(R));
  w = 2 .* (multiconv(R, D) - X)';
  for k = 1:K
    cost_grad(:, k) = w * sptoeplitz([D(:, k); zeros(size(R, 1) - 1,1)], [D(1, k); zeros(size(R, 1) - 1,1)]);
  end

  R_grad = cost_grad;

end
