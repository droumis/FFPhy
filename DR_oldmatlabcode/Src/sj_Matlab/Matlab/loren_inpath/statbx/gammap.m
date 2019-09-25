function p = gammap(q,a)
%GAMMAP Gamma distribution function.
%	GAMMAP(Q,A) is Pr(X < Q) where X is a Gamma random variable with
%	shape parameter A.  A must be scalar.
%
%	See also GAMMAQ, GAMMAR.

%  GKS  2 August 1998.

if a < 0, error('Gamma parameter A must be positive'); end
p = gammainc(q,a);
k = find(q < 0);
if any(k), p(k) = zeros(size(k)); end
