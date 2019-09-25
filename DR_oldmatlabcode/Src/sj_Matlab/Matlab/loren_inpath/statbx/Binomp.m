function f = binomp(x,n,p);
%BINOMP	BINOMP(X,N,P) is the cumulative distribution function at X of the
%	binomial distribution with denominator N and success probability P.
%	X, N and P must be scalars.

% GKS 13 Aug 95, 6 Aug 98

if n <= 0
   error('Binomial denominator must be positive.')
elseif (p < 0) | (p > 1)
   error('Binomial probability must be >= 0 and <= 1.')
elseif x >= n
   f = 1;
elseif x < 0
   f = 0;
else
   x = floor(x);
   f = betainc( 1-p, n-x, x+1 );
end
