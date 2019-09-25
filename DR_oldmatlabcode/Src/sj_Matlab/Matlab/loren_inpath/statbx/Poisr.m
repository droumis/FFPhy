function y = poisr(lambda,n)
%POISR Poisson random number generator.
%	POISR(LAMBDA,N) generates N random deviates from the
%	Poisson distribution with mean LAMBDA.
%
%	POISR(LAMBDA) generates a single random deviate.
%
%	See POISP, RAND.

% GKS 31 July 1993, 28 July 1999
% Uses naive inversion method.

if nargin==1, n=1; end;
low = max( 0 , floor( lambda - 8*sqrt(lambda) ));
hi = ceil( lambda + 8*sqrt(lambda) + 4/lambda );
x = (low:hi)';
p = cumsum( exp( -lambda + x.*log(lambda) - gammaln(x+1) ));
u = rand(n,1);
y = zeros(n,1);
for i = 1:n,
   k = find( u(i) < p );
   y(i) = x( k(1) );
end;
