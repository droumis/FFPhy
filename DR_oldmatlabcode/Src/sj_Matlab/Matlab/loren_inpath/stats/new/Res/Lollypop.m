% LOLLYPOP: Produces 3D scatterplot in which points are anchored onto the 
%           X,Y plane by a stem (subtended line).  Does not open a new figure window.
%           Calls lollipop().
%
%     Syntax:  lollypop(x,y,z,{color})  OR  lollypop([x,y,z],{color})
%
%           x,y,z = column vectors of identical length.
%           color = optional character value indicating color or circle and line
%                     [default = 'k'].
%

% RE Strauss, 2/18/98
%   6/6/03 - added optional color argument.

function lollypop(x,y,z,color)
  if (nargin < 3)
    if (nargin < 2)
      color = [];
    else
      color = y;
    end;
    z = x(:,3);
    y = x(:,2);      
    x = x(:,1);
  elseif (nargin < 4)
    color = [];
  end;
  
  lollipop(x,y,z,color);
  return;
  