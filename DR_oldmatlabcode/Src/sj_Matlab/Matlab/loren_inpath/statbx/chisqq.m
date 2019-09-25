function q = chisqq(p,v)
%CHISQQ	Quantiles of the chi-square distribution.
%	Q = CHISQQ(P,V) satisfies P(X < Q) = P, where X follows a
%	chi-squared distribution on V degrees of freedom.
%	V must be a scalar.
%
%	See also CHISQP.

%	Gordon K Smyth, University of Queensland, gks@maths.uq.edu.au
%	27 July 1999

%	Reference:  Johnson and Kotz (1970). Continuous Univariate
%	Distributions, Volume I. Wiley, New York.

q = 2*gammaq(p,v/2);
