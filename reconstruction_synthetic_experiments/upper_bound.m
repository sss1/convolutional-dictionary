function ub = upper_bound(N, n, lambda, sigma)

  ub = 4 .* lambda .* sigma .* sqrt(2 .* n .* log(2 .* N)) ./ N;

end
