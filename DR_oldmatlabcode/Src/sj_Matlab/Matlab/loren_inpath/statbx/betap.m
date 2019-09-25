function p = betap(x,a,b);
%BETAP	Beta cumulative distribution function.
%	BETAP(X,A,B) gives P(X < x) where X follows a BETA distribution
%	with parameters A > 0 and B > 0.  A and B must be scalars.
%
%	See also BETAQ, BETAR

%	Gordon Smyth, gks@maths.uq.edu.au, University of Queensland
%	31 July 1998

p = zeros(size(x));

k = find(x >= 1);
if any(k), p(k) = ones(size(k)); end

k = find(x > 0 & x < 1);
if any(k), p(k) = betainc(x(k),a,b); end
   
