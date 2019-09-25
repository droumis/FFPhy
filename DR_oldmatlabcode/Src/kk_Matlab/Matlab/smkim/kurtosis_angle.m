function k = kurtosis_angle(a,dim)
%KURTOSIS_ANGLE Circular kurtosis of sample(s).
%
%   K = KURTOSIS_ANGLE(A) takes angles (in radians) in A and computes the sample
%   circular kurtosis. For vector input, K is a scalar. For matrix input, K is a
%   row vector containing the circular kurtosis of each column of A. For N-D
%   arrays, KURTOSIS_ANGLE operates along the first non-singleton dimension.
%
%   KURTOSIS_ANGLE(A,DIM) takes the circular kurtosis along the dimension DIM of
%   A.
%
%Written by SMK, 2009 December 9
%

if (exist('nansum') ~= 2)
  error(['DISPERSION_ANGLE depends on the m-file NANSUM ' ...
      '(MATLAB Statistics Toolbox)']);
end
if (exist('nanmean') ~= 2)
  error(['DISPERSION_ANGLE depends on the m-file NANMEAN ' ...
      '(MATLAB Statistics Toolbox)']);
end

if ~isfloat(a) || ~isreal(a) || ~all(isfinite(a(:)))
  error('A must be a floating-point real array with Inf elements');
end
if any((a(:) < -pi) | (a(:) > +pi))
  %warning('Some values in A lie outside the interval [-pi,+pi]');
end

if (nargin < 2) || isempty(dim)
  % The output size for [] is a special case, handle it here.
  if isequal(a,[])
    d = nan(class(a)); 
    return; 
  end;
  % Figure out which dimension to work along
  dim = find(size(a) ~= 1,1,'first');
  if isempty(dim)
    dim = 1;
  end
end
n = size(a,dim);

% Take mean angle along specified dimension, excluding NaN values
z = nanmean(complex(cos(a),sin(a)),dim);
Rbar = abs(z);
a_mean = angle(z);
a_centered = bsxfun(@minus_angle,a,a_mean);

% Real component of second central moment
alpha2 = real(nanmean(complex(cos(2*a_centered),sin(2*a_centered)),dim));

k = (alpha2 - Rbar.^4) ./ (1 - Rbar).^2;

