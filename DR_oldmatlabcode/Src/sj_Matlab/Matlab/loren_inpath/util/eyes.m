function a = eyes(d, n, m)
% eyes - array of block identity matrices
%
%	EYES(D) is the same as EYE(D).  EYES(D, N, M) returns a matrix
%	composed of N by M copies of the D-order identity matrix.
%
%	See also: EYE, ONES, ZEROS

usage = char([9 'EYES(D, N, M)']);

if (nargin == 0) disp(usage); return; end

if (nargin == 3)
  a = zeros(n*d, m*d);
  a(1:d, 1:d) = eye(d);
  for r = 0:n-1
    a(r*d + (1:d), 1:d) = a(1:d, 1:d);
  end
  for c = 1:m-1
    a(:, c*d + (1:d)) = a(:, 1:d);
  end
else
  a = eye(d);
end