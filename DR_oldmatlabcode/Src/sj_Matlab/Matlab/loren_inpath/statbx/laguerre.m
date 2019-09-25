function [t,w]=laguerre(n,alpha)
%LAGUERRE Calculate nodes and weights for generalized laguerre
%	quadrature.  [T W]=LAGUERRE(N,ALPHA) puts N nodes in T and
%	corresponding weights in W so that the integral from 0 to
%	infty of f(x) x^alpha exp(-x) is approximated by W' * f(T).
%	See also LAGUERR2.

%	Adapted from Netlib routine gaussq.f
%	GKS  15 May 92

i = (1:n)';
a = 2 .* i - 1 + alpha;
i = (1:(n-1))';
b = sqrt( i .* (i + alpha) );

[v d] = eig( diag(a) + diag(b,1) + diag(b,-1) );
w = v(1,:)';
[t i] = sort( diag(d) );
w = w(i);

w = gamma(alpha+1) .* w.^2;
