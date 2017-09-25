function ub = upper_bound(N, n, lambda, sigma, p)

  ub = 4 .* lambda .* sigma .* N.^((1 - p)/p) .* n.^(max(0,(p - 2)/(2*p)));

end
