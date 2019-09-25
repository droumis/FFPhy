% PROJECT:  Given a set of points specified by vectors x,y and a line 
%           specified by its slope and intercept, projects the points onto 
%           the line, returning a single column vector of centered scores.
%
%       Usage: scores = project(x,y,slope,intcpt,{nocenter})
%
%         x,y =       vectors of point coordinates for n points.
%         slope =     slope of reference line.
%         intcpt =    intercept of reference line.
%         nocenter =  optional boolean flag indicating, if true, that scores
%                       are not to be centered.
%         ------------------------------------------------------------------
%         scores =   [n x 1] vector of scores.
%

% RE Strauss, 1/20/98
%   3/14/02 - added 'nocenter' option, improved documentation.

function scores = project(x,y,slope,intcpt,nocenter)
  if (nargin < 5) nocenter = []; end;
  
  if (isempty(nocenter))
    nocenter = 0;
  end;
  
  lenx = length(x);
  leny = length(y);
  if (lenx ~= leny)
    error('  Project: vectors of X,Y coordinates must be same length');
  end;
  scores = zeros(lenx,1);

  l1 = -1/slope;
  l2 = -1;
  l3 = 0;
  dl = sqrt(l1*l1 + l2*l2);

  scores = -(x*l1 + y*l2 + l3)/dl;
  if (~nocenter)
    scores = scores - mean(scores);
  end;

  return;
