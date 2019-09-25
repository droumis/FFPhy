function r = chisqr(m,n,v)
%CHISQR	Random deviates from the chi-square distribution.
%	CHISQR(M,N,V) generates an MxN matrix of random deviates from
%	the chi-squared distribution on V degrees of freedom.
%	V must be a scalar.
%
%	See also CHISQP, CHISQQ.

%	Gordon K Smyth, University of Queensland, gks@maths.uq.edu.au
%	9 Dec 1999

%	Reference:  Johnson and Kotz (1970). Continuous Univariate
%	Distributions, Volume I. Wiley, New York.

r = zeros(m,n);
for i=1:m,
for j=1:n,
   r(i,j) = 2*gammar(v/2);
end
end
