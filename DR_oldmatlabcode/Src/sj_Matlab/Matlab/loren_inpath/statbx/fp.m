function p = fp(q,v1,v2)
%FP	F cumulative distribution function.
%	FP(Q,V1,V2) gives Pr(X < Q) where X follows an F distribution
%	on V1 and V2 degrees of freedom.
%
%	See also FQ, FR.

%	Gordon Smyth, University of Queensland, gks@maths.uq.edu.au
%	5 Jun 92, 28 Jan 93, 3 Apr 97, 1 Aug 98.

%	Reference:  Johnson & Kotz, Chapter 26, page 78.

a = v1 ./ 2;
b = v2 ./ 2;
q = v1 .* q ./ ( v2 + v1 .* q) ;
p = betap(q,a,b);
