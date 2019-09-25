function q = tq(p,v);
%TQ	t distribution quantiles.
%	Q = TQ(P,V) satisfies Pr(T < Q) = P where T follows a
%	t-distribution on V degrees of freedom.
%	V must be a scalar.

%	Gordon Smyth, University of Queensland, gks@maths.uq.edu.au
%	2 August 1998

if v <= 0, error('Degrees of freedom must be positive.'); end;

q = sign(p-0.5)*sqrt( fq(2*abs(p-0.5),1,v) );
