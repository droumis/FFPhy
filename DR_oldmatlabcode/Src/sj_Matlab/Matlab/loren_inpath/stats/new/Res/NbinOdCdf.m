% NBINODCDF: Cdf for the negative binomial distribution, using the
%             mean-variance (m,v) parameterization.
%
%     Usage: Y = nbinodcdf(X,m,v)
%
%         X = integer values of X at which the function is to be evaluated.
%         m = mean parameter.
%         v = variance parameter.
%         -----------------------------------------------------------------
%         Y = distribution values evaluated at X.
%

% Hilborn, R. & M. Mangel. 1997. The ecological detective: confronting
%   models with data.  Princeton Univ. Press.

% RE Strauss, 2/19/03

function Y = nbinodcdf(X,m,v)
  if (nargin==0) help nbinodcdf; return; end;
  if (nargin<3)
    error('  NBINODCDF: too few input arguments.');
  end;
  
  if (any(~isintegr(X(:))))
    error('  NBINODCDF: X must consist of integer values.');
  end;
  
  Y = zeros(size(X));
  [r,c] = size(X);
  
  k = (m*m)/(v-m);
  
  for ir = 1:r
    for ic = 1:c
      x = X(ir,ic);
      t1 = gamma(k+x)/(gamma(k)*prod(1:x));
      t2 = (k/(k+m))^k;
      t3 = (m/(m+k))^x;
      Y(ir,ic) = t1*t2*t3;
    end;
  end;
  
  return;
  
