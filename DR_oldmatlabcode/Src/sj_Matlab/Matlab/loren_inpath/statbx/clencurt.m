function [Q,error] = clencurt(funfcn,a,b,n)
%CLENCURT Numerical evaluation of an integral, variable point method.
%	Q = CLENCURT('f',a,b) approximates the integral of f(x) from a to b
%	using an 11-point Clenshaw-Curtis formula.
%
%	Q = CLENCURT('f',a,b,n) uses n+1-point Clenshaw-Curtis integration.
%	[Q,error] = CLENCURT('f',a,b,n) computes a (usually conservative)
%	estimate of the approximation error.
%	
%	The integral is exact for polynomials of degree n or less.

%	GKS  23 Jan 1997

if nargin<4,
   n = 10;
else
   % make sure n is even
   n = 2*ceil(n/2);
end;
s = (0:n)';
s2 = 2*(0:(n/2))';
x = cos(pi*s/n);
f = feval(funfcn,x*(b-a)/2+(b+a)/2);
f(1) = f(1)/2;
f(n+1) = f(n+1)/2;
c = 2/n * cos(s2*s'*pi/n) * f;
c(1) = c(1)/2;
c(n/2+1) = c(n/2+1)/2;
Q = (b-a)*sum( -c./(s2-1)./(s2+1) );
error = 2*max(abs(c(n/2)),abs(c(n/2+1)));
