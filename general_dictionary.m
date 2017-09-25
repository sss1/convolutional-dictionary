function [R, D, opt] = learn_constrained_sparse_dictionary(X, n, K, varargin)
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
%   lambda_R_sparse (real scalar > 0,     Default: 0.01) weight of sparsity penalty on R
%   gamma_0         (real scalar >= 0,    Default: 0.01) initial gradient descent step size
%   num_iterations  (integer scalar >= 0, Default: 100)  total number of gradient descent iterations
% 
% Outputs:
%   R (non-negative, real, sparse, (N - n + 1) X K): matrix whose jth column encodes positions at which feature j occurs
%   D (real, n X K): matrix whose j^th column encodes the j^th feature
%   opt (non-negative, real, num_iterations X 1): history of objective function

  % default values for optional arguments
  lambda_R_sparse = 0.01;
  gamma_0 = 0.01;
  num_iterations = 200;

  % parse optional arguments
  for i=1:2:length(varargin)
     switch varargin{i}
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
     otherwise
       error(sprintf('%s is not a valid argument name', varargin{i}));
     end
  end

  [N, d] = size(X);

  % Initialize each D_j as a uniformly random unit vector (with 0's at each end)
  D = zeros(n, K, d);
  for dim_idx = 1:d
    D(:, :, dim_idx) = normc(normrnd(0, 1, [n K]));
  end

  % Initialize each R_j as a non-negative IID gaussian vector
  R = max(normrnd(0, var(X)./K, [(N - n + 1), K]), 0);

  if nargout >= 3
    opt = zeros(num_iterations, 1);
  end

  for iteration = 1:num_iterations
    gamma = gamma_0 .* (iteration.^(-1/2));

    % Update D
    D = D - gamma * grad_D(D, R, X);
    for dim_idx = 1:d
      D(:, :, dim_idx) = normc(D(:, :, dim_idx)); % ensure that columns of D are unit vectors
    end

    % Update R
    R = R - gamma * grad_R(D, R, X);
    R = max(R, 0); % ensure that R is non-negative
    R(:) = ProjectOntoL1Ball(R(:), lambda_R_sparse); % ensure that ||R||_1 <= lambda_R_sparse

    % Compute objective function
    if nargout >= 3
      for i = 1:d
        opt(iteration) = opt(iteration) + norm(multiconv(R, D(:, :, i)) - X(:, i));
      end
      opt(iteration) = opt(iteration) + norm(R(:), 1);
    end

  end

end

% Helper function: calculates the gradient of the objective with respect to D
function D_grad = grad_D(D, R, X, lambda_D_smooth, zero_edge)

  [n, K, d] = size(D);

  % Gradient of reconstruction error
  D_grad = zeros(size(D));

  w = zeros(size(X));
  for dim_idx = 1:d
    w(:, dim_idx) = multiconv(R, D(:, :, dim_idx)) - X(:, dim_idx);
  end
  for k = 1:K
    % For a column vector c, sptoeplitz([c; zeros(n - 1,1)], zeros(n,1))
    % is a sparse equivalent of convmtx(c, n).
    R_convmtx = spconvmtx(R(:, k), n);
    D_grad(:, k, :) = w' * R_convmtx;
  end

end

% Helper function: calculates the gradient of the objective with respect to R
function R_grad = grad_R(D, R, X)

  [n, K, d] = size(D);

  % Gradient of reconstruction error
  R_grad = zeros(size(R));
  for dim_idx = 1:d
    w = multiconv(R, D(:, :, dim_idx)) - X(:, dim_idx);
    for k = 1:K
      D_convmtx = spconvmtx(D(:, k, dim_idx), size(R, 1));
      R_grad(:, k) = R_grad(:, k) + D_convmtx' * w;
    end
  end

end

function C = spconvmtx(v, n)
  C = sptoeplitz([v; zeros(n - 1,1)], [v(1); zeros(n - 1,1)]);
end
