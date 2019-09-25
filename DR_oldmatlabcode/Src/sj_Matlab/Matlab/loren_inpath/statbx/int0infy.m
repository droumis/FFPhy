function int=int0infy(fun,n,alpha)
%INT0INFY Integrate function from 0 to infty by laguerre quadrature.
%	INT0INFY('FUN',N,ALPHA) integrates FUN(x) by approximating
%	FUN(x) x^(-ALPHA) EXP(x) with an Nth order polynomial.
%	ALPHA is assumed zero if omitted and N defaults to 32 if omitted.

%	GKS  24 Dec 92

if nargin<2, n=32; end;
if nargin<3, alpha=0; end;
[t w] = laguerre(n,alpha);
int = w'* ( feval(fun,t) .* t.^(-alpha) .* exp(t) );
