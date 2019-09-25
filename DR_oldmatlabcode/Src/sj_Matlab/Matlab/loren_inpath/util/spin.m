function frames = spin(origin, delta, ax)
% spin - generate animations of rotation.
%
% 	M = SPIN() returns a movie describing a spin of the current
%	plot objects in 10 degree increments about a vertical axis through
%	[0,0].  The movie can be played using MOVIE.  
%
%	M = SPIN(ORIGIN, DELTA, AXIS) spins the objects about the
%	specified 3-vector origin and axis (specified as [AZ, EL]) in
%	increments of DELTA.  Each of the arguments is optional.  See
%	ROTATE.
%
%	SPIN rotates from the current view and in the current axis.
%	These should therefore be set correctly using VIEW and AXIS
%	before calling spin.


if (nargin < 1) origin = [0,0,0]; end;

if (nargin < 2) delta  = 10;	  end;

if (nargin < 3) ax     = [0,90];  end;

% Freeze axes
axis(axis);

nframes = 360/delta;
frames = moviein(nframes);

frames(:, 1) = getframe;

handles = get(gca, 'Children');

for i = 2:nframes
  for h = handles'
    rotate(h, ax, delta, origin);
  end
  frames(:,i) = getframe;
end

