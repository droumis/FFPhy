% GEOGDIST: Geographic (great-circle arc) distance among a set of sites
%           specified by latitude and longitude.  If latitude and/or
%           longitude are specified by vectors of up to length 3, the elements
%           are assumed to represent degrees, minutes, and seconds.
%           Northern latitudes are positive, southern latitudes are negative.
%           Eastern longitudes are positive, western longitudes are negative.
%           Only the degrees portion of the coordinate need be specified as
%           negative; signs of minutes and seconds are ignored.
%
%     Usage: dist = geogdist(lat,long,{units})
%
%           lat =   [n x 1] (degrees) or [n x 2] (degrees, minutes) or 
%                     [n x 3] (degrees, minutes, seconds) matrix of latitudinal 
%                     coordinates for n locations.
%           long =  matching matrix of longitudinal coordinates for n points.
%           units = optional units of output distances [default = 1 = km]:
%                     1 = distance in kilometers
%                     2 = distance in miles
%                     3 = distance in nautical miles (=minutes)
%           -------------------------------------------------------------------
%           dist =  [n x n] matrix of pairwise geographic distances.
%

% Ayres Jr., F.  1954.  Theory and Problems of Plane and Spherical Trigonometry.
%   Schaum's Outline Series in Mathematics.  McGraw-Hill.  Chapter 23.

% RE Strauss, 6/29/96

function dist = geogdist(lat,long,units)
  if (nargin < 3) units = []; end;

  if (isempty(units))
    units = 1;                        % Default units = kilometers
  end;

  [n1,p1] = size(lat);
  [n2,p2] = size(long);

  if (n1 ~= n2)
    error('GEOGDIST: numbers of locations must be identical for latitude, longitude');
  else
    n = n1;
  end;

  if (p1 > 1)                         % Convert latitude minutes, seconds to hundredths
    s = sign(lat(:,1));               %   of degree
    c = abs(lat(:,1)) + abs(lat(:,2)./60);
    if (p1 > 2)
      c = c + abs(lat(:,3)./3600);
    end;
    lat = s.*c;
  end;

  if (p2 > 1)                         % Convert longitude minutes, seconds to hundredths
    s = sign(long(:,1));              %   of degree
    c = abs(long(:,1)) + abs(long(:,2)./60);
    if (p2 > 2)
      c = c + abs(long(:,3)./3600);
    end;
    long = s.*c;
  end;

  dist = zeros(n,n);                  % Allocate output distance matrix

  for i = 1:(n-1)                     % Cycle thru all possible pairs of localities
    lat1 = lat(i);
    long1 = long(i);

    for j = (i+1):n
      lat2 = lat(j);
      long2 = long(j);

      if (mean([lat1 lat2]) < 0)          % Reflect southern to northern hemisphere
        lat1 = -lat1;
        lat2 = -lat2;
      end;

      if (mean([long1 long2]) < 0)        % Reflect western to eastern hemisphere
        long1 = -long1;
        long2 = -long2;
      end;

      a = 2*pi*(90-lat2)/360;             % Convert colatitudes to radians
      b = 2*pi*(90-lat1)/360;

      if (a==b)                         % Same latitude?
        C = abs(long1-long2)*60;          % Convert longitude diff to minutes
        nmi = abs(C * cos((pi/2)-a));     % Convert to minutes (=nautical miles)
      else                              % Different latitudes?
        C = 2*pi*abs(long1-long2)/360;    % Convert longitude diff to radians
        AB1 = atan(cos((b-a)/2)*sec((b+a)/2)*(1/(tan(C/2)+eps)));
        AB2 = atan(sin((b-a)/2)*(1/sin((b+a)/2))*(1/(tan(C/2)+eps)));
        A = AB1 - AB2;
        B = 2*AB1 - A;
        c = 2 * atan(tan((b-a)/2)*sin((B+A)/2)*(1/(sin((B-A)/2)+eps)));
        c = c*360/(2*pi);                  % Convert radians to degrees
        nmi = abs(c*60);                   % Convert to minutes (=nautical miles)
      end;

      dist(i,j) = nmi;                     % Stash distance as nautical miles
      dist(j,i) = nmi;
    end;
  end;

  if (units == 1)
    dist = dist * 1.852;                   % Convert to kilometers
  elseif (units == 2)
    dist = dist * 1.151;                   % Convert to statute miles
  end;

  return;
