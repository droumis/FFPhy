% NormalTest: Measures consistency with a normal distribution based 
%             on the Shapiro-Francia W' test, a modification of the Shapiro-Wilks' 
%             test (Royston 1993).  The test statistic W' is approximately the 
%             squared correlation between the sorted data and the corresponding 
%             rankits or normal scores.
%               Handles censored or uncensored data; censored data are passed as 
%             NaN values. 
%               The randomized significance level of the statistic is based on the 
%             null hypothesis that the sample has been drawn from a normally 
%             distributed population; thus the test is left-tailed.  
%             For censored data, the proportion of censored values is held constant
%             during randomization.
%             
%     Usage: [W,prob,a,xa] = normaltest(x,{iter},{censdir})
%
%         x =         data vector.
%         iter =      optional number of iterations for randomized significance level
%                       [default = 5000].
%         censdir =   optional boolean variable indicating the direction of censoring:
%                       -1 = left-censored [default];
%                       +1 = right-censored.
%         ----------------------------------------------------------------------------
%         W =       test statistic value.
%         prob =    significance level of W under the null hypothesis of a
%                     normally distributed population.
%         a =       Shapiro-Francia coefficients.
%         xa =      sorted non-missing data values corresponding to coefficients.
%

% Royston, P. 1993. A toolkit for testing for non-normality in complete and
%   censored samples.  The Statistician 42:37-43.

% RE Strauss, 12/10/02
%   1/2/03 - output coefficients and matching data.

function [W,prob,a,xa] = normaltest(x,iter,censdir)
  if (nargin < 1) help normaltest; return; end;

  if (nargin < 2) iter = []; end;
  if (nargin < 3) censdir = []; end;
  
  if (nargout < 2) iter = 0; end;
  
  if (isempty(iter))    iter = 5000; end;
  if (isempty(censdir)) censdir = -1; end;
  
  ntot = length(x);
  i = find(isfinite(x));
  ngood = length(i);
  ncens = ntot-ngood;
  x = x(i);

  if (ncens>0 & censdir<0)          % Convert left-censored values to right-censored
    x = -x;
  end;
  x = sort(x(:));                   % Sort input data
  xa = x;

  if (ngood < 4)
    error('  NormalTest: non-censored N must be > 3.');
  end;
  
  m = rankits(ntot);                % Mean order statistics
  a = (m'*m)^(-0.5)*m;              % Coefficients
  a = a(1:ngood);
  W = corr(a,x).^2;
  
  prob = NaN;
  if (iter)                         % Randomized probability
    prob = 0;
    incr = 1/iter;
    for it = 1:iter
      xr = sort(randn(ntot,1));
      xr = xr(1:ngood);
      Wr = corr(a,xr).^2;
      if (Wr <= W)
        prob = prob + incr;
      end;
    end;
  end;
  
  if (ncens>0 & censdir)            % Reverse coefficients if necessary
    a = -a;
  end;
  
  return;
