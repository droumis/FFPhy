function r = fr(m,n,v1,v2)
%FR	F distribution random deviates.
%	FR(M,N,V1,V2) generates an M by N matrix of random deviates
%	from the F distribution on V1 and V2 degrees of freedom.
%	V1 and V2 must be scalars.
%
%	See also FP, FQ

%	Gordon Smyth, gks@maths.uq.edu.au, University of Queensland

r = betar(m,n,v1/2,v2/2);
r = v2/v1 * r./(1-r);
