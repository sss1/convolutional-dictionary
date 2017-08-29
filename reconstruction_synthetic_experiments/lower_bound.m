function lb = lower_bound(N, n, lambda, sigma)


  t1 = lambda/(8 .* N);
  t2 = lambda .* n;
  t3 = sigma .* sqrt(log(N - n + 1));

  lb = t1 .* min(t2, t3);

end
