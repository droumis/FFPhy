% CIRCFIT:  Fits a circle to a set of 2-dimensional points by finding the center point
%           having minimum variance in Euclidean distance to all points.
%           Fits exact circles to 2 and 3 (noncolinear) points.
%
%     Usage: [center,radius] = circfit(crds,{doplot})
%
%           crds =    [n x 2] matrix of point coordinates.
%           doplot =  optional boolean flag indicating whether (=1) or not (=0) the 
%                       points and circle are to be plotted [default = 0].
%           -------------------------------------------------------------------------
%           center =  [1 x 2] row vector of point coordinates of fitted center.
%           radius =  radius of circle.
%

% RE Strauss, 2/24/98
%   9/7/99 -  changed plot colors for Matlab v5, plus misc changes.
%   10/7/99 - changed call to circcrds().
%   3/21/00 - changed call to circcrds().

function [center, radius] = circfit(crds,doplot)
  if (nargin < 2) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  [n,p] = size(crds);

  if (n==2 & p~=2)                % Transpose input matrix if necessary
    crds = crds';
    [n,c] = size(crds);
  end;

  if (n<2)
    error('  CIRCFIT: need minimum of 2 points to fit a circle');
  end;
  if (p~=2)
    error('  CIRCFIT: 2-dimensional points only');
  end;

  if (n==3)
    if (isonline(crds(1,:),[crds(2,:),crds(3,:)],1))
      error('  CIRCFIT: 3 points must not be colinear');
    end;
  end;

  center = mean(crds);            % Use centroid as initial estimate of center

  center = fmins('circfitf',center,[],[],crds);
  radius = mean(eucl(crds,center));

  if (doplot)
    figure;
    plot(crds(:,1),crds(:,2),'ok');
    ccrds = circcrds(radius,center);
    hold on;
    plot(ccrds(:,1),ccrds(:,2),'b');
    hold off;
    sqplot([crds;ccrds]);
  end;

  return;
