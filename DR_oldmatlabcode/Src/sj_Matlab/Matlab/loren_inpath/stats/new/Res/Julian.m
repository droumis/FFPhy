% JULIAN: Given the year, month and day (and, optionally, hours, minutes and 
%         seconds), returns the Julian date (day of the year, 1-365 or 1-366).  
%         Values greater than the normal range are carried to the next unit 
%         (e.g., month values greater than 12 are carried to years).
%         Two-digit year is considered to be in the current century.  Handles 
%         all leap years.
%
%     Usage: jdate = julian(y,m,d,{h},{mi},{s},{series})
%                   OR
%            jdate = julian(date,{series})
%
%         y,m,d,h,mi,s = values for year, month, day, hour, min, sec.  These 
%                          may be vectors of identical length k; if a mixture 
%                          of vectors and scalars, the scalars are expanded.
%         date =         vector of [y,m,d,{h},{mi},{s}], or matrix with k such 
%                          rows.
%         series =       optional boolean flag indicating, if true, that the 
%                          input dates represent a sequential series, and that 
%                          Julian dates are to be accumulated across years.
%         ----------------------------------------------------------------------
%         jdate =        vector of julian date(s).
%

% RE Strauss, 2/23/99.  Modified from Matlab function datenum().
%   3/13/99 -   added series option.
%   12/30/00 -  miscellaneous improvements.

function jdate = julian(y,m,d,h,mi,s,series)
  if (nargin < 4) h = []; end;
  if (nargin < 5) mi = []; end;
  if (nargin < 6) s = []; end;
  if (nargin < 7) series = []; end;

  if (nargin < 3)                     % Expand single matrix into col vectors
    if (nargin == 2)                    % 2-argument form of input
      series = m;
    end;

    [y,m,d,h,mi,s] = extrcols(y);
  end;

  if (isempty(m) | isempty(d))
    error('  JULIAN: too few input arguments');
  end;

  if (isempty(series))
    series = 0;
  end;

  if (isempty(h))  h = 0;  end;
  if (isempty(mi)) mi = 0; end;
  if (isempty(s))  s = 0;  end;

  if (isscalar(y))  y =  y(ones(size(s)));  end;  % Expand scalars into vectors
  if (isscalar(m))  m =  m(ones(size(y)));  end;
  if (isscalar(d))  d =  d(ones(size(m)));  end;
  if (isscalar(h))  h =  h(ones(size(d)));  end;
  if (isscalar(mi)) mi = mi(ones(size(h))); end;
  if (isscalar(s))  s =  s(ones(size(mi))); end;
  if (isscalar(y))  y =  y(ones(size(s)));  end;

  i = find(y<100);
  if (~isempty(i))                    % Expand 2-digit year
    i = find(y<100 & y>20);             % Cutoff for century 1900
    if (~isempty(i))
      y(i) = 1900 + y(i);
    end;

    i = find(y<=20);                    % Cutoff for century 2000
    if (~isempty(i))
      y(i) = 2000 + y(i);
    end;
  end;

  i = find(m==0);                     % Check range of month
  if (~isempty(i))
    m(i) = ones(length(i),1);
  end;
  i = find(m>12);
  while (~isempty(i))
    m(i) = m(i) - 12;
    y(i) = y(i) + 1;
    i = find(m>12);
  end;
                                      % Running total of days per month
  cumdpm = cumsum([0;31;28;31;30;31;30;31;31;30;31;30;31]); 

  % result = (365 days/year)*(number of years) + number of leap years
  %      + days in previous months + days in this month + fraction of a day

  ndate = length(y);
  jdate = zeros(ndate,1);

  for i = 1:ndate
    jdate(i) = 365.*y(i) ...          % Convert year, month, day to date number
      + ceil(y(i)/4)-ceil(y(i)/100)+ceil(y(i)/400) ...
      + reshape(cumdpm(m(i)),size(m(i))) ...
      + ((m(i)>2) & ((rem(y(i),4)==0 & rem(y(i),100)~=0) | rem(y(i),400)==0)) ...
      + d(i);

    jdate(i) = jdate(i) + (h(i).*3600+mi(i).*60+s(i))./86400;

    y(i) = y(i)-1;                    % Standardize to date of current year
    m(i) = 12;
    d(i) = 31;

    zd = 365.*y(i) ...
      + ceil(y(i)/4)-ceil(y(i)/100)+ceil(y(i)/400) ...
      + reshape(cumdpm(m(i)),size(m(i))) ...
      + ((m(i)>2) & ((rem(y(i),4)==0 & rem(y(i),100)~=0) | rem(y(i),400)==0)) ...
      + d(i);

    jdate(i) = jdate(i) - zd;
  end;

  if (series)
    pos = find(jdate(1:ndate-1)>jdate(2:ndate))+1;
    displ = 0;
    while(~isempty(pos))
      p1 = pos(1);
      if (length(pos)>1)
        p2 = pos(2)-1;
      else
        p2 = ndate;
      end;
      displ = displ + julian([y(p1-1)+1 12 31]);
      jdate(p1:p2) = jdate(p1:p2) + displ;
      pos = find(jdate(1:ndate-1)>jdate(2:ndate))+1;
    end;
  end;

  return;
