function d = dispersion_angle(a,dim)
%DISPERSION_ANGLE Sample circular dispersion of angle data
%
%   D = DISPERSION_ANGLE(A) takes angles (in radians) in A and computes the
%   sample circular dispersion. For vector input, D is a scalar. For matrix
%   input, D is a row vector containing the circular dispersion of each column
%   of A. For N-D arrays, DISPERSION_ANGLE operates along the first
%   non-singleton dimension.
%
%   DISPERSION_ANGLE(A,DIM) takes the circular dispersion along the dimension
%   DIM of A.
%
%   DISPERSION_ANGLE treats NaNs as missing values, and removes them.
%
%Depends on:
%   MINUS_ANGLE (written by SMK)
%   NANSUM (MATLAB Statistics Toolbox)
%   NANMEAN (MATLAB Statistics Toolbox)
%
%Written by SMK, 2009 December 9
%
%TODO: Is this the correct way to compute circular dispersion? According to
%Fisher (1993), a von Mises distribution with concentration parameter kappa
%should have a circular dispersion of
%
% d = 1 ./ kappa .* besseli(0,kappa) ./ besseli(1,kappa);
%
%But applying DISPERSION_ANGLE to simulated random von Mises samples does not
%produce results that agree with this theoretical value.
%

warning('This function is not fully tested and may not give correct results.');

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
z = nansum(complex(cos(a),sin(a)),dim);
a_mean = angle(z);

% Reshape and tile to align with A to allow subtraction
repsz = ones([1 ndims(a)]);
repsz(dim) = size(a,dim);
a_centered = minus_angle(a,repmat(a_mean,repsz));

% Length of first central moment
rho1 = abs(nanmean(complex(cos(a_centered),sin(a_centered)),dim));

% Real component of second central moment
alpha2 = real(nanmean(complex(cos(2*a_centered),sin(2*a_centered)),dim));

d = 0.5*(1 - alpha2) ./ (rho1.^2);


