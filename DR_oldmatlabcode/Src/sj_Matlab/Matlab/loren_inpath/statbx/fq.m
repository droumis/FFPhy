function q = fq(p,v1,v2)
%FQ	F distribution quantiles.
%	Q = FQ(P,V1,V2) satisfies Pr(X < Q) = P, where X follows
%	and F distribution on V1 and V2 degrees of freedom.
%	V1 and V2 must be scalars.
%
%	See also FP

%	Gordon Smyth, gks@maths.uq.edu.au, University of Queensland

q = betaq(p,v1/2,v2/2);
q = v2/v1 * q./(1-q);
