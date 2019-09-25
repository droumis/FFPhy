function c = minus_angle(a,b)
%MINUS_ANGLE Difference between angles, accounting for circular wraparound
%
%   C = MINUS_ANGLE(A,B) computes the angular difference A - B, where A and B
%   are interpreted as containing angles (in radians). A and B must have the
%   same dimensions unless one is a scalar. The return value C is wrapped to be
%   bounded between -pi and +pi, and can be interpreted as the smallest angle by
%   which B must be rotated to match A.
%
%Written by SMK, 2009 October 15.
%

if ~isfloat(a) || ~isreal(a)
  error('A must be floating-point real');
end
if any((a < -pi) | (a > pi))
  %warning('some A values are outside the interval [-pi,+pi]');
end

if ~isfloat(b) || ~isreal(b)
  error('B must be floating-point real');
end
if any((b < -pi) | (b > pi))
  %warning('some A values are outside the interval [-pi,+pi]');
end

try
  c = angle(exp(i*(a - b)));
catch
  error('Non-scalar arguments must match in size.');
end


