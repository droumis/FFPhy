% PROBITF: Objective function for probit(); returns the negative of the 
%           log-likelihood function, given the data and a vector or parameter 
%           values.
%
%     Usage: negL = probitf(b,x,y)
%
%         b = row vector of estimates of coefficients: [b1 b0].
%         x = [n x p] matrix of independent variables.
%         y = [n x 1] vector of binary dependent variable.
%         ----------------------------------------------------------------------
%         negL = negative of the log-likelihood function.  If -L is minimized, then 
%               L is maximized.
%

% RE Strauss, 4/2/99

function negL = probitf(b,x,y)
  p = normcdf(x*b');                          % Predicted probabilities

  tol = 1e-9;
  if (min(p)<tol)
    i = find(p<tol);
    p(i) = tol*ones(length(i),1);
  end;
  if (max(p)>(1-tol))
    i = find(p>(1-tol));
    p(i) = (1-tol)*ones(length(i),1);
  end;

  negL = -sum(y.*log(p) + (1-y).*log(1-p));   % Neg log-likelihood
%b_minp_maxp_negL = [b min(p) max(p) negL]

  return;

