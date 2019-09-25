function c = sum_angle(a,b)
%SUM_ANGLE Sum of angles, accounting for circular wraparound
%
%   C = SUM_ANGLE(A,B) computes the angular sum A + B, where A and B are
%   interpreted as containing angles (in radians). A and B must have the same
%   dimensions unless one is a scalar. The return value C is wrapped to be
%   bounded between -pi and +pi, and can be interpreted as the resultant angle
%   of rotation after rotating by A and B.
%
%Written by SMK, 2009 October 15.
%

if ~isfloat(a) || ~isreal(a)
  error('A must be floating-point real');
end
if any((a < -pi) | (a > pi))
  warning('some A values are outside the interval [-pi,+pi]');
end

if ~isfloat(b) || ~isreal(b)
  error('B must be floating-point real');
end
if any((b < -pi) | (b > pi))
  warning('some A values are outside the interval [-pi,+pi]');
end

try
  c = angle( complex(cos(a),sin(a)) .* complex(cos(b),sin(b)) );
catch
  error('Non-scalar arguments must match in size.');
end

