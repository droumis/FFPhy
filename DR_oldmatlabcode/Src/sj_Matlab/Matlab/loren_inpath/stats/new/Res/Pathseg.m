% PATHSEG:  Given a set of coordinates describing a path in p dimensions, and a given 
%           number of straight line segments, finds the line-segment length required 
%           to exactly step off the line segments along the path.  Returns the 
%           line-segment length and the coordinates on the path corresponding to the 
%           steps.  
%
%           Note that there may be >1 optimum to this problem.  This function works by 
%           starting with a segment length shorter than required, and increasing the 
%           segment length until the stepping process overshoots the path.  Repeats 
%           this process, systematically decreasing the setment-length increment, 
%           until a minimum residual (not necessarily zero) is attained.
%
%           Removes any points with non-finite coordinate values from path, 
%           effectively connecting the points on either side by a straight line.
%
%     Usage: [steplen,stepcrds,resid] = pathseg(path,nsteps,{doplot})
%
%           path =     [n x p] coordinates of the path, for n points in p dimensions.
%           nsteps =   number of steps of the line segment along the path.
%           doplot =   optional boolean flag indicating, if true, that a plot of the 
%                        path and step coordinates is to be produced [default = 0].
%           --------------------------------------------------------------------------
%           steplen =  length of line segment required to exactly step off the path.
%           stepcrds = coordinates which the line segments intersect the path.
%           resid =    residual distance left over at the end of the path.
%

% RE Strauss, 4/17/98
%   6/1/01 -  removes points with missing (non-finite) coordinate values.
%   6/29/01 - corrected problem with unitialized variables.

function [steplen,stepcrds,resid] = pathseg(path,nsteps,doplot)
  if (nargin < 3) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  i = find(~isfinite(rowsum(path)));          % Remove non-finite points
  if (~isempty(i))
    path(i,:) = [];
  end;

  [n,p] = size(path);
  tol = 1e-7;

  d = path(2:n,:) - path(1:n-1,:);            % Length of input path
  pathlen = sum(sqrt(sum((d.*d)')));

  Lsteplen = 0.5 *( pathlen./nsteps);         % Initial guess at lower bound of step length 
  delta = 2 * Lsteplen;
  min_resid = pathlen;
  min_resid_steplen = Lsteplen;

  while (delta > tol)
    delta = 0.1 * delta;
    resid = -0.1;
    steplen = Lsteplen;

    while (resid < 0)
      resid = pathsegf(path,steplen,nsteps);
      if (resid <= 0)
        Lsteplen = steplen;
      end;
      if (abs(resid) < abs(min_resid))
        min_resid = resid;
        min_resid_steplen = steplen;
      end;
      steplen = steplen + delta;
    end;
  end;

  steplen = min_resid_steplen;

  [resid,stepcrds] = pathsegf(path,steplen,nsteps,doplot);

  return;

