% MedianCI: Finds an exact confidence interval (>= a specified alpha-level) 
%           for the median, by column.
%
%     Usage: [m,ci,rcl] = medianci(X,alpha)
%
%         X =     [r x c] data matrix.
%         alpha = optional confidence level [default = 95].
%         ---------------------------------------------------------------------
%         m =     [1 x c] vector of medians, for c columns.
%         ci =    [2 x c] corresponding matrix of critical values for the medians.
%         rcl =   realized confidence level (min >= alpha).
%

% RE Strauss, 5/31/02

function [m,ci,rcl] = medianci(X,alpha)
  if (nargin < 2) alpha = []; end;
  
  if (isvector(X))
    X = X(:);
  end;
  [n,c] = size(X);
  
  if (isempty(alpha))
    alpha = 0.95;
  end;
  if (alpha > 1)
    alpha = alpha/100;
  end;
  
  p = (1-alpha)/2;                          % 1-tailed probability
  z = norminv(1-p);                         % Normal approx to notch depths
  d = floor(0.5*(n+1 - z*sqrt(n)));         % Notch depth
  d = max([d,1]);
  rcl = 1-2*binocdf(d-1,n,0.5);             % Realized confidence level
  
  while (rcl > alpha)
    d = d+1;
    rcl = 1-2*binocdf(d-1,n,0.5);
  end;
  
  m = zeros(1,c);                           % Allocate output matrices
  ci = zeros(2,c);
  
  for ic = 1:c                              % Cycle thru columns
    x = sort(X(:,ic));                        % Isolate and sort column
    md = (n+1)/2;
    if (isintegr(md))                         % Median
      m(ic) = x(md);
    else
      m(ic) = 0.5*(x(floor(md))+x(ceil(md)));
    end;
    ci(:,ic) = [x(d) x(n-d+1)]';
  end;
  
  

  return;
  