% HISTARW:  Returns position of arrow and adjusted axis ranges, for histogram 
%           data.
%
%     Usage: [position,shaftlen,v] = histarw(x,n,centers,v,shaftlen)
%
%           x =         [n x 1] vector of data values, or scalar value of horizontal 
%                         position of tip of arrow (x-value of arrow).
%           n =         vector of counts per bin.
%           centers =   corresponding centers of bins.
%           v =         axis ranges: [xmin xmax ymin ymax].
%           shaftlen =  relative (<1) or absolute (>1) height of arrow shaft.
%           -----------------------------------------------------------------
%           position =  4-element vector of coordinates of arrow:
%                         [xtip,ytip,xbase,ybase].
%           shaftlen =  absolute height of arrow shaft.
%           v =         adjusted axis ranges.
%           

% RE Strauss, 11/25/99
%   2/26/00 - fixed problem with calculation of ymax.
%   2/ 8/01 - pass either data or value of mean.

function [position,shaftlen,v] = histarw(x,n,centers,v,shaftlen)
  if (isscalar(x))                          % Mean of data
    xtip = x;
  else
    xtip = mean(x);                           
  end;
  [cm,i] = min(abs(centers-xtip));          % Find bin containing mean

  binheight = n(i);                         % Height of bar in this bin
  if (i-1 > 0)
    binheight = max([binheight, n(i-1)]);
  end;
  if (i+1 <= length(n))
    binheight = max([binheight, n(i+1)]);
  end;

  ymax = v(4);                              % Current y-max

  if (shaftlen < 1)                         % Convert relative shaft length
    shaftlen = ymax*shaftlen;               %   to absolute
  end;

  ytip = binheight + 0.25*shaftlen;         % Ordinate of arrow tip
  ybase = ytip + shaftlen;                  % Ordinate of arrow base

  ymax = max([ymax, ybase+(0.05*ybase)]);   % New y-maximum
  v(4) = ymax;

  position = [xtip ytip xtip ybase];        % Arrow position

  return;
