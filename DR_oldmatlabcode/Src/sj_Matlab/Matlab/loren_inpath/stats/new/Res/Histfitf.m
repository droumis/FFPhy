% Histfitf: Objective function for histfit. Finds sum-of-squared deviation of 
%           predicted from observed values, subject to the constraint that 
%           predicted values must be non-negative.
%
%     Usage: sse = histfitf(b,[d1 d2 d3 d4],freqs)
%

% RE Strauss, 2/25/01

function sse = histfitf(b,devs,freqs)
  n = size(devs,1);
  pred = [ones(n,1) devs]*b;

  e = (pred-freqs).^2;
  i = find(pred<0);
  e(i) = 10*e(i);

  sse = sum(e);

  return;