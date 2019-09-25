% GETPTID: Given an existing scatterplot (in the active graphics window) and the 
%          coordinates of the points that have been plotted, allows the user to 
%          identify particular points by digitizing their positions on the plot.
%
%     Usage:
%           ptids = getptid(crds)
%                     or
%           ptids = getptid(x,y)
%
%           crds =  [n x 2] matrix of point coordinates.
%           x,y, =  matching vectors (length n) of point coordinates.
%           ---------------------------------------------------------
%           ptids = indices of digitized points.
%

% RE Strauss, 3/24/98
%   6/8/00 - accomodate change in Matlab v5 ginput function; misc improvements.

function ptids = getptid(x,y)
  if (nargin==1)
    [n,p] = size(x);
    if (p~=2)
      error('  GETPTID: input matrix must be 2-dimensional');
    end;
    y = x(:,2);
    x = x(:,1);
  end;

  disp(' ');
  disp('  Begin digitizing points.');
  disp('  Press <Enter> key (with cursor in window) to end input.');

  ptids = [];                           % Allocate return matrix
  goodpoint = 1;
  while (goodpoint)                       % Capture points
    [xx,yy,button] = ginput(1);
    if (isempty(button))                  % End of points
      goodpoint = 0;
    end;

    if (goodpoint)
      dist = eucl([x y],[xx yy]);           % Distances of digitized pt to plotted pts
      [mindist,i] = min(dist);              % Find minimum distance
      ptids = [ptids; i];
    end;
  end;

  return;

