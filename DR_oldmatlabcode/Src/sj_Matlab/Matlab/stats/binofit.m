function [phat, pci] = binofit(x,n,alpha)
%BINOFIT Parameter estimates and confidence intervals for binomial data.
%   PHAT = BINOFIT(X,N) Returns estimates of the probability of success for
%   the binomial distribution.  X and N are scalars containing the number of
%   successes and the number of trials, respectively.  If X and N are vectors,
%   BINOFIT returns a vector of estimates whose I-th element is the parameter
%   estimate for X(I) and N(I).  A scalar value for X or N is expanded to the
%   same size as the other input.
%
%   [PHAT, PCI] = BINOFIT(X,N,ALPHA) gives MLEs and 100(1-ALPHA) 
%   percent confidence intervals given the data. Each row of PCI contains
%   the lower and upper bounds for the corresponding element of PHAT.
%   By default, the optional parameter ALPHA = 0.05 corresponding to 95%
%   confidence intervals.
%
%   Note:  Unlike most other distribution fitting functions, BINOFIT treats a
%   vector X as a collection of measurements from separate samples, and
%   returns a vector of estimates.  If you want to treat X as a single sample
%   and compute a single parameter estimate and confidence interval, use
%   BINOFIT(SUM(X),SUM(N)) when N is a vector, and BINOFIT(SUM(X),N*LENGTH(X))
%   when N is a scalar.
%
%   See also BINOCDF, BINOINV, BINOPDF, BINORND, BINOSTAT, MLE. 

%   Reference:
%      [1]  Johnson, Norman L., Kotz, Samuel, & Kemp, Adrienne W.,
%      "Univariate Discrete Distributions, Second Edition", Wiley
%      1992 p. 124-130.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 2.16.2.6 $  $Date: 2006/04/25 07:19:41 $

if nargin < 3 
    alpha = 0.05;
end


% Initialize params to zero.
[row, col] = size(x);
if min(row,col) ~= 1
   error('stats:binofit:VectorRequired','First argument must be a vector.');
end

[r1,c1] = size(n);
if ~isscalar(n)
   if row ~= r1 || col ~= c1
      error('stats:binofit:InputSizeMismatch',...
            'The first two inputs must match in size.');
   end
end
if ~isfloat(x)
   x = double(x);
end

if any(n<0) || any(n~=round(n)) || any(isinf(n)) || any(x>n)
    error('stats:binofit:InvalidN',...
          'All N values must be non-negative integers at least as large as X.')
end
if any(x<0)
    error('stats:binofit:InvalidX','X must not be negative.')
end
phat = x./n;

if nargout > 1,
   nu1 = 2*x;
   nu2 = 2*(n-x+1);

   F   = finv(alpha/2,nu1,nu2);
   lb  = (nu1.*F)./(nu2 + nu1.*F);

   xeq0 = find(x == 0);
   if ~isempty(xeq0)
      lb(xeq0) = 0;
   end
   
   nu1 = 2*(x+1);
   nu2 = 2*(n-x);
    
   F   = finv(1-alpha/2,nu1,nu2);
   ub = (nu1.*F)./(nu2 + nu1.*F);
   
   xeqn = find(x == n);
   if ~isempty(xeqn)
      ub(xeqn) = 1;
   end

   pci = [lb(:) ub(:)];
end

