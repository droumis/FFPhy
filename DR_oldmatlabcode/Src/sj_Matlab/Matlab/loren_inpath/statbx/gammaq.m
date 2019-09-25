function q = gammaq(p,a);
%GAMMAQ	Gamma distribution quantiles.
%	Q = GAMMAQ(P,A) satisfies Pr(X < Q) = P where X follows a
%	Gamma distribution with shape parameter A > 0.
%	A must be scalar.
%
%	See also GAMMAP, GAMMAR

%	Gordon Smyth, gks@maths.uq.edu.au, University of Queensland
%	2 August 1998

%	Method:  Newton iteration on the scale of log(X), starting
%	from point of inflexion of cdf.  Monotonic convergence is
%	guaranteed.

if any(p < 0 | p > 1), error('P must be between 0 and 1'); end

if a == 0.5, q = 0.5*normq((1-p)/2).^2; return; end
if a == 1, q = -log(1-p); return; end

G = gammaln(a);
x = log(a);
q = a;
for i=1:10,
    x = x - (gammap(q,a) - p) ./ exp(a*x - q - G);
    q = exp(x);
end

k = find((p>0 & p<1.5e-4) | (p>0.99999 & p<1));
if any(k),
   X = x(k);
   Q = q(k);
   P = p(k);
   for i=1:10,
      X = X - (gammap(Q,a) - P) ./ exp(a*X - Q - G);
      Q = exp(X);
   end;
   q(k) = Q;
end

k = find(p == 0);
if any(k), q(k) = zeros(size(k)); end

k = find(p == Inf);
if any(k), q(k) = ones(size(k)); end
