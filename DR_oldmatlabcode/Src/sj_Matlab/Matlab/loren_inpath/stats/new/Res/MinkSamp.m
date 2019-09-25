% MINKSAMP: Sample from a Minkowski distribution.  Provide as input either a 
%           sample size for a random sample, or a vector of values on the 
%           interval [0,1] for a systematic sample.
%
%     Usage: x = minksamp(n,k)
%
%           n = sample size, or vector of values on interval [0,1].
%           k = Minkowski coefficient.
%           -------------------------------------------------------
%           x = sample from Minkowski distribution.
%

% RE Strauss & EG Dyreson, 5/24/01

function y = minksamp(n,k)
  n = n(:);

  if (isscalar(n))
    u = rand(n,1);
  else
    u = n;
  end;

  n = 1000;

  x = linspace(1,15,n)';
  x = x-mean(x);
  y = exp(-(abs(x).^k)/2);   
  
  cdf = cumsum(y)./sum(y);              % Cumulative distribution

  y = zeros(size(u));
  for i = 1:length(u)
    d = cdf - u(i);
    [m,indx] = max(d(d<0));

    x1 = x(indx);                       % Linearly interpolate
    x2 = x(indx+1);
    y1 = d(indx);
    y2 = d(indx+1);

    y(i) = x2 - y2*(x2-x1)/(y2-y1);
  end;

  return;
