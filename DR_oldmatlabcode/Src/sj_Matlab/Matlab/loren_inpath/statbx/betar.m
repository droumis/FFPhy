function r = betar(m,n,a,b);
%BETAR	Beta distribution random deviates.
%	BETAR(M,N,A,B) generates an M x N matrix of random deviates
%	from the BETA distribution with parameters A > 0 and B > 0.
%	A and B must be scalars.
%
%	See also BETAP, BETAQ

%	Gordon Smyth, gks@maths.uq.edu.au, University of Queensland
%	1 August 1998

r = betaq(rand(m,n),a,b);
