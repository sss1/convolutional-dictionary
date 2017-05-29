% Given two matrices A and B with equal width, the jth column of C is the
% convolution of the jth column of A with the jth column of B
function C = multiconv(A, B)
  [Na, Wa] = size(A);
  [Nb, Wb] = size(B);
  if Wa ~= Wb
    error(['A and B must have the same width. ' ...
           'Currently, A has width ' num2str(Wa) ...
           ' and B has width ' num2str(Wb) '.']);
  end
  C = zeros(Na + Nb - 1, Wa);
  for col_idx = 1:Wa
    C(:, col_idx) = conv(A(:, col_idx), B(:, col_idx));
  end
end
