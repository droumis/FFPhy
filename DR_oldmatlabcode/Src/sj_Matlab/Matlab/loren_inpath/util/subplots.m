function hout = subplots (m, n, p)
% SUBPLOTS -- create subplots
% 	SUBPLOTS(M, N) calls SUBPLOT M*N times in order to create all the
% 	subplots needed to form an M by N grid.
%
%	SUBPLOTS(M, N, [A B ...]) creates only the subplots specified
%	by A, B etc.  Subplots are numbered row first, beginning with
%	1 in the upper left hand corner.
%
%	H = SUBPLOTS(...) returns a vector of object handles pointing
%	to the axes that are created.
%
%	See also SUBPLOT, AXES.

if (nargin < 3) p = 1:m*n; end

h = [];

for pp = p
  h = [h, subplot(m,n,pp)];
end

if (nargout > 0) hout = h; end
