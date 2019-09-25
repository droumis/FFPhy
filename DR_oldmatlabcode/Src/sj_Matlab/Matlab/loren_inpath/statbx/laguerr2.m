function [t,w]=laguerr2(n,alpha)
%LAGUERR2 [T W]=LAGUERR2(N,ALPHA) calculates N nodes T and weights W
%	in order to approximate the expectation of a function of a
%	gamma random variable by generalized laguerre quadrature.
%	If X has a standard gamma distribution with shape parameter
%	ALPHA, then E[ g(X) ] is approximated by W' * g(T).
%	If X has scale parameter BETA, then E[ g(X) ] is approximated
%	by W' * g(BETA * T).  See LAGUERRE.

%	GKS  15 May 92

i = (1:n)';
a = 2 .* i - 2 + alpha;
i = (1:(n-1))';
b = sqrt( i .* (i + alpha - 1) );

[v d] = eig( diag(a) + diag(b,1) + diag(b,-1) );
w = v(1,:)';
[t i] = sort( diag(d) );
w = w(i).^2;
