function p = tp(t,v);
%TP	TP(T,V) is the cumulative distribution function at T of the
%	t-distribution on V degrees of freedom.
%	V must be a scalar.

% Gordon Smyth, University of Queensland, gks@maths.uq.edu.au
% 3 Apr 97

if v <= 0, error('Degrees of freedom must be positive.'); end;
x = t.^2 ./ (v + t.^2) ;
p = 0.5 .* ( 1 + sign(t) .* betainc( x, 0.5, 0.5*v ) ) ;
