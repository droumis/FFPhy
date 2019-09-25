% BETANC: Returns cumulative probability of x (0<x<1) for the noncentral beta
%         distribution with parameters a,b and noncentrality lambda.
%
%     Syntax: prob = betanc(x,a,b,lambda)
%
%         x =      scalar or matrix of abscissa limits to be evaluated
%         a,b =    scalars or corresponding matrices of shape parameters of
%                    beta distribution
%         lambda = scalar or corresponding matrix of noncentrality parameter(s)
%         ===================================================================
%         prob =   matrix of probabilities corresponding to x
%

% Lenth, R.V. 1987.  Algorithm AS 226: Computing noncentral beta probabilities.
%   Appl. Stat. 36:241-244.
% Frick, H. 1990. Algorithm AS R84: A remark on Algorithm AS 226: Computing
%   noncentral beta probabilities.  Appl. Stat. 39:311-312.

% Implementation notes:
%   Frick notes that, for extremely large values of lambda, the number of
% iterations can be very large.  This was found to be the case when applied to
% confidence intervals for Mahalanobis distances.  For this reason Length's
% recommended maxiter=200 was increased to 1000 (max observed during tests was
% 440).  However, Frick's adjustment of the starting index m was found to
% produce cumulative probabilities > 1 when the index adjustment kicked in, and
% so was abandoned.

function prob = betanc(x,a,b,lambda)

  % Check dimensions of input arguments

  if (nargin < 4)
    lambda = 0;
  end;

  [rows,cols] = size(x);
  matsize = rows * cols;
  x = x(:);

  [r,c] = size(a);
  if (r*c ~= matsize)
    if (matsize>1 & r*c==1)
      a = a * ones(matsize,1);
    else
      error('  Error: input matrices incompatible');
    end;
  else
    a = a(:);
  end;

  [r,c] = size(b);
  if (r*c ~= matsize)
    if (matsize>1 & r*c==1)
      b = b * ones(matsize,1);
    else
      error('  Error: input matrices incompatible');
    end;
  else
    b = b(:);
  end;

  [r,c] = size(lambda);
  if (r*c ~= matsize)
    if (matsize>1 & r*c==1)
      lambda = lambda * ones(matsize,1);
    else
      error('  Error: input matrices incompatible');
    end;
  else
    lambda = lambda(:);
  end;

  if (any(x<0 | x>1))
    x
    error('  Error: x out of range [0,1]');
  end;

  if (any(a<=0 | b<=0 | lambda<0))
    error('  Error: parameters out of range');
  end;

  % Estimate cumulative probabilities

  maxiter = 1000;
  maxerr = 1e-6;

  prob = zeros(matsize,1);

  for i = 1:matsize
    if (x(i)<maxerr)              % Check extremes of range of x
      prob(i) = 0;
    elseif (x(i)>(1-maxerr))
      prob(i) = 1;
    else                          % Else find prob of non-extreme value

      xi = x(i);
      ai = a(i);
      bi = b(i);
      lambdai = lambda(i);

      % Initialize the series

      c = lambdai/2;
      beta = gammaln(ai) + gammaln(bi) - gammaln(ai+bi);
      temp = betainc(xi,ai,bi);
      gx = exp(ai*log(xi) + bi*log(1-xi) - beta - log(ai));
      q = exp(-c);

      it = 0;
      errbound = maxerr+1;

      ax = q * temp;
      sumq = 1 - q;
      prob(i) = ax;

      % Iterate over subsequent terms

      while ((it<maxiter) & (errbound>maxerr))
        it = it + 1;
        temp = temp - gx;
        gx = xi * (ai+bi+it - 1) * gx / (ai+it);
        q = q*c/it;
        sumq = sumq - q;
        ax = temp * q;
        prob(i) = prob(i) + ax;
        errbound = (temp-gx) * sumq;
      end;
    end;
  end;

  return;

