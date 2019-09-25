function q = betaq(p,a,b);
%BETAQ	Beta distribution quantiles.
%	Q = BETAQ(P,A,B) satisfies Pr(X < Q) = P where X follows a
%	BETA distribution with parameters A > 0 and B > 0.
%	A and B must be scalars.
%
%	See also BETAP

%	Gordon Smyth, smyth@wehi.edu.au,
%	Walter and Eliza Hall Institute of Medical Research
%	1 August 1998

%	Method:  Work with the distribution function of log(X/(1-X)).
%	The cdf of this distribution has one point of inflexion at the
%	the mode of the distribution at log(A/B).  Newton's method
%	converges monotonically from this starting value.

q = zeros(size(p));

k = find(p >= 1);
if any(k), q(k) = ones(size(k)); end

k = find(p > 0 & p < 1);
if any(k),
   P = p(k);
   B = betaln(a,b);
   r = a/(a+b);
   x = log(a/b);
   ex = a/b;
   for i=1:6,
      x = x - (betap(r,a,b) - P) ./ exp(a*x - (a+b)*log(1+ex) - B);
      ex = exp(x);
      r = ex./(1+ex);
   end
   q(k) = r;
end
