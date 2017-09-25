function ub = upper_bound(N, n, lambda, sigma, is_independent)

  ub = 4 .* lambda .* sigma .* sqrt(2 .* n .* log(2 .* N)) ./ N;
  if nargin >= 5 && is_independent % when noise is independent, use the tighter bound for this case
    ub = min(ub, 8 .* lambda .* sigma .* sqrt(2 .* log(2 .* N)) ./ N);
  end

end
