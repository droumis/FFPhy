% PTSCALE: Scales the coordinates of a point on the current plot 
%          to or from proportions of the distance along the abscissa 
%          and ordinate, based on the axis minima and maxima.
%
%     Syntax: spt = ptscale(x,y,dir)
%
%          pt =  2-element vector of point coordinates.
%          dir = direction of the mapping [default=0]:
%                   0 = from true coords to proportional coords.
%                   1 = from proportional coords to true coords.
%          -------------------------------------------------------
%          spt = corresponding vector of scaled point coordinates.
%

% RE Strauss, 2/7/97

function spt = ptscale(pt,dir)
  if (nargin < 2) dir = []; end;

  if (isempty(dir)) dir = 0; end;

  spt = pt;
  bnds = axis;                  % Boundaries of current plot axes

  if (~dir)                     % Transform true to proportional
    spt(1) = (pt(1)-bnds(1)) / (bnds(2)-bnds(1));
    spt(2) = (pt(2)-bnds(3)) / (bnds(4)-bnds(3));
  end;
  
  if (dir)                      % Transform proportional to true
    spt(1) = bnds(1) + pt(1)*(bnds(2)-bnds(1));
    spt(2) = bnds(3) + pt(2)*(bnds(4)-bnds(3));
  end;

  return;

