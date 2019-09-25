% POLARCRD: Transforms polar coordinates to rectangular coordinates, 
%           or vice versa.
%
%     Usage: [r,theta] = polarcrd(x,y,0)  [default]
%                     OR
%            [x,y] = polarcrd(r,theta,1)
%
%            r,theta, x,y = vectors of identical length (except that
%                             one of input vectors may be a scalar).
%

% RE Strauss, 9/16/97
%   3/19/00 - allow single-matrix input: [x, y];
%             corrected errors for x==0 and y==0.

function [a,b] = polarcrd(c,d,polar_to_rect)
  if (nargin < 2) d = []; end;
  if (nargin < 3) polar_to_rect = []; end;
  
  if (isempty(d))
    if (size(c,2)>1)
      d = c(:,2);
      c = c(:,1);
    else
      error('  POLARCRD: too few input arguments.');
    end;
  end;

  if (isempty(polar_to_rect))
    polar_to_rect = 0;
  end;

  lengthc = length(c);
  lengthd = length(d);

  if (lengthc==1 & lengthd>1)
    c = c * ones(size(d));
  end;
  if (lengthd==1 & lengthc>1)
    d = d * ones(size(c));
  end;

  if (polar_to_rect)              % Polar to cartesian
    a = c .* cos(d);                % a = x
    b = c .* sin(d);                % b = y
  else                            % Cartesian to polar
    a = sqrt(c.*c + d.*d);          % a = r
    i = find(a<eps);
    if (~isempty(i))
      a(i) = eps*ones(length(i));
    end;

    b = asin(d./a);                 % b = theta
    if (~isempty(i))
      a(i) = zeros(length(i));
    end;

    for i = 1:length(b)             % Adjust for quadrant
      if (c(i)<0 & d(i)>0)
        b(i) = pi - b(i);
      elseif (c(i)<0 & abs(d(i))<eps)
        b(i) = pi;
      elseif (c(i)<0 & d(i)<0)
        b(i) = abs(b(i)) + pi;
      elseif (c(i)>0 & d(i)<0)
        b(i) = b(i) + 2*pi;
      elseif (abs(c(i))<eps & d(i)<0)
        b(i) = 3*pi/2;
      end;
    end;
  end;

  return;
