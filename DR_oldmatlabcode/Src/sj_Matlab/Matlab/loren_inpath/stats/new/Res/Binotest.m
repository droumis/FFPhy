% BINOTEST: Take replicate samples from a known binomial distribution with N=3 
% (simulating HW proportions), and % find (1) p by enumeration, 
% and (2) p from the best-fitting (least-squares) binomial distribution.

N = 2;
n = 10;

iter = 10;
p1_iter = zeros(iter,1);
p2_iter = zeros(iter,1);
results = [];

X = [0:N]';

for true_p = 0.10:0.10:0.90
  for it = 1:iter
    b = binornd(N,true_p,n,1);      % Random column vector
    [freq,x] = hist(b,X);
    p1_iter(it) = sum(freq'.*X)/(n*N);
    p2_iter(it) = fmin('binofitf',0,1,[],N,freq'); % Optimize p
  end;

  p1 = mean(p1_iter);
  p1_std = std(p1_iter);

  p2 = mean(p2_iter);
  p2_std = std(p2_iter);

[n true_p p1 p1_std p2 p2_std]
  results = [results; n true_p p1 p1_std p2 p2_std];
end;

results

