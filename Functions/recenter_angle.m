function b = recenter_angle(a,center)
%RECENTER_ANGLE Remap angles to an arbitrary 2*pi interval.
%
%   B = RECENTER_ANGLE(A,CENTER) takes a floating-point array of angles in
%   radians and maps the values mod (2*pi) to the interval [CENTER-pi,
%   CENTER+pi]. B approximately preserves the trigonometric properties of A,
%   subject to floating-point error on the order of
%   max(max(eps(A)),max(eps(B))).
%
%   RECENTER_ANGLE(A) is the same as RECENTER_ANGLE(A,0).
%
%Written by SMK, 2009 October 15.
%

if (nargin == 1)
  center = 0;
end

if ~isfloat(a) || ~isreal(a)
  error('A must be floating-point real');
end
if any((a < -pi) | (a > pi))
  %warning('some A values are outside the interval [-pi,+pi]');
end

if ~isscalar(center) || ~isfloat(center) || ~isreal(center) || ...
    ~isfinite(center)
  error('CENTER must be floating-point finite real scalar');
end

% It's very important to perform all computations in double precision, because
% MATLAB's mod function accumulates floating-point errors
center = double(center);
b = mod(double(a) + pi - center,2*pi) - pi + center;

b = cast(b,class(a));

