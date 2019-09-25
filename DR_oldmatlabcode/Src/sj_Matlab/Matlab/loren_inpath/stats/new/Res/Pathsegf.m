% PATHSEGF: Given a set of coordinates describing a path in p dimensions, a line-segment 
%           length, and a specified number of straight line segments, steps off the line 
%           segments along the path and returns the residual.  The residual is negative 
%           if the stepping process stops short of the end of the path, positive if it 
%           over-steps the path.  In the latter case, the last few step coordinates are 
%           repeated.
%
%     Usage: [resid,stepcrds] = pathsegf(path,steplen,nsteps,{doplot})
%
%           path =     [n x p] coordinates of the path, for n points in p dimensions.
%           steplen =  length of line segment.
%           nsteps =   number of steps of the line segment along the path.
%           doplot =   optional boolean flag indicating, if true, that a plot of the 
%                        path and step coordinates is to be produced [default = 0].
%           --------------------------------------------------------------------------
%           resid =    residual distance left over at the end of the path.
%           stepcrds = coordinates which the line segments intersect the path.
%

% RE Strauss, 4/11/98
%   8/21/99 - changed plot colors for Matlab v5.

function [resid,stepcrds] = pathsegf(path,steplen,nsteps,doplot)
  if (nargin < 4) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;
  
  [n,p] = size(path);
  tol = 1e-5;

  savepath = path;

  stepcrds = zeros(nsteps+1,p);               % Initialize output crd matrix
  stepcrds = path(1,:);

  for step = 1:nsteps
    pathdist = path - ones(n,1)*stepcrds(step,:);     % Dists of path crds from current step
    pathdist = sqrt(sum((pathdist.^2)')') - steplen;  % Resid dists from current step

    if (any(pathdist>0))
      i = min(find(pathdist>0));                      % Closest path crd with positive resid

      A = stepcrds(step,:);
      B = path(i-1,:);
      C = path(i,:);
      c = sqrt(sum((B-C).*(B-C)));
      pl = 0;
      pu = c;
  
      if (A==B)                                         % If next intsect on next path segment,
        Pl = A;
        Pu = C;
        P = (Pl+Pu)./2;
        p = (pl+pu)./2;
        a = sqrt(sum((P-A).*(P-A))) - steplen;
  
        while (abs(a) > tol)                              % Binary search for intersection point
          if (a>0)
            Pu = P;
            pu = p;
          else
            Pl = P;
            pl = p;
          end;
          P = (Pl+Pu)./2;
          p = (pl+pu)./2;
          a = sqrt(sum((P-A).*(P-A))) - steplen;
        end;
    
      else                                              % If next intsect not on next path segment,
        cos_alpha = cos(pi - abs(angledev(A,B,C)));
        b = pathdist(i-1) + steplen;
    
        Pl = B;
        Pu = C;
        P = (Pl+Pu)./2;
        p = (pl+pu)./2;
        
        a = sqrt(b*b + p*p - 2*b*p*cos_alpha) - steplen;  % Use law of cosines to find dist to path
        while (abs(a) > tol)                              % Binary search for intersection point
          if (a>0)
            Pu = P;
            pu = p;
          else
            Pl = P;
            pl = p;
          end;
          P = (Pl+Pu)./2;
          p = (pl+pu)./2;
          a = sqrt(b*b + p*p - 2*b*p*cos_alpha) - steplen;
        end;
      end;
    
      stepcrds(step+1,:) = P;
      path(1:i-1,:) = ones(i-1,1)*stepcrds(step+1,:);

    else
      stepcrds(step+1,:) = stepcrds(step,:);
    end;
  end;
  
  d = stepcrds(nsteps+1,:)-path(n,:);       % Residual path length
  resid = sqrt(sum(d.*d));                  % Pos if undershoot, neg if overshoot
  if (any(pathdist>0))
    resid = -resid;
  end;

  if (doplot)
    figure;
    plot(savepath(:,1),savepath(:,2),'k');
    hold on;
    plot(stepcrds(:,1),stepcrds(:,2),'ko');
    hold off;
    putbnd(savepath(:,1),savepath(:,2));
    axis('equal');
  end;

  return;

