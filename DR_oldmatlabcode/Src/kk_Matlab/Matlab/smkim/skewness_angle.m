function s = skewness_angle(a,dim)
%SKEWNESS_ANGLE Circular skewness of sample(s).
%
%   S = SKEWNESS_ANGLE(A) takes angles (in radians) in A and computes the sample
%   circular skewness. For vector input, S is a scalar. For matrix input, S is a
%   row vector containing the circular skewness of each column of A. For N-D
%   arrays, SKEWNESS_ANGLE operates along the first non-singleton dimension.
%
%   SKEWNESS_ANGLE(A,DIM) takes the circular skewness along the dimension DIM of
%   A.
%
%Depends on:
%   MINUS_ANGLE (written by SMK)
%   NANSUM (MATLAB Statistics Toolbox)
%   NANMEAN (MATLAB Statistics Toolbox)
%
%Written by SMK, 2009 December 9
%

if (exist('minus_angle') ~= 2)
  error(['DISPERSION_ANGLE depends on the m-file MINUS_ANGLE ' ...
      '(written by SMK)']);
end
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
a_mean = angle(z);
Rbar = abs(z);
a_centered = bsxfun(@minus_angle,a,a_mean);

% Angle of first central moment
mu1 = angle(nanmean(complex(cos(a_centered),sin(a_centered)),dim));

% Imaginary component of second central moment
beta2 = imag(nanmean(complex(cos(2*a_centered),sin(2*a_centered)),dim));

s = beta2 ./ power(sqrt(1 - Rbar),3);




