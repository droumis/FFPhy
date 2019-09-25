function p = chisqp(q,v)
%CHISQP Chi-square distribution function.
%	CHISQP(Q,V) is the probability that a chi-squared random
%	variable on V degrees of freedom is < Q.  V must be a scalar.
%
%	See also CHISQQ.

%  GKS  5 June 1992.  Revised 2 August 1998.

p = gammainc(q/2, v/2);
