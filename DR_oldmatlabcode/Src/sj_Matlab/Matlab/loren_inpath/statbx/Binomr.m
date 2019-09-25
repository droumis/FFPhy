function y = binomr(n,p)
%BINOMR Binomial random number generator.
%	BINOMR(N,P) generates a single random deviate from the binomial
%	distribution with denominator n and success probability P.
%
%	See RAND.

% GKS 12 Aug 1995, 28 Jul 1999
% Uses naive inversion method.

f = zeros(n+1,1);
for i=0:n, f(i+1) = binomp(i,n,p); end
u = rand(1);
k = find( u < f );
y = k(1) - 1;
