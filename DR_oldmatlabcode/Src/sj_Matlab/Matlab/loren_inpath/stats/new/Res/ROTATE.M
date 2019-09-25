% ROTATE: Given a configuration of points, rotates a 2D configuration about 
%         a given point by a specified angle (in radians).
%
%     Usage: [newpts,pivot] = rotate(pts,theta,{pivot},{doplot})
%
%         pts =    [N x 2] matrix of coordinates of point configuration.
%         theta =  angle by which configuration is to be rotated; a positive 
%                    angle rotates the configuration counterclockwise, a 
%                    negative angle rotates it clockwise.
%         pivot =  optional 2-element vector of coordinates of the pivot point
%                    [default = centroid].
%         doplot = optional boolean variable indicating, if true, that plots
%                    are to be produced depicting the point configuration before
%                    and after rotation [default = 0].
%         ----------------------------------------------------------------------
%         newpts = [N x 2] matrix of registered & rotated points.
%         pivot =  [1 x 2] vector of coordinates of pivot point.
%

% RE Strauss, 6/26/96
%   5/6/03 - return pivot point.
%   5/14/03 - added optional plots.

function [newpts,pivot] = rotate(pts,theta,pivot,doplot)
  if (~nargin) help rotate; return; end;
  
  if (nargin < 3) pivot = []; end;
  if (nargin < 4) doplot = []; end;
  
  if (isempty(doplot)) doplot = 0; end;

  N = size(pts,1);                        % Number of points

  if (isempty(pivot))
    [area,perim,pivot] = polyarea(pts);   % Use centroid for pivot
  else
    pivot = pivot(:)';
    if (length(pivot)~=2)
      error('  Rotate: pivot point must be vector of length 2.');
    end;
  end;  

  savepts = pts;
  pts = pts - ones(N,1)*pivot;            % Zero-center on pivot
  dev = anglerotation([0 0],[1,0],pts,1); % Angular deviations of pts from horizontal
  dev = dev + theta;                      % Add angle of rotation to deviations
  r = sqrt(pts(:,1).^2 + pts(:,2).^2);    % Distances of pts from origin

  newpts = zeros(size(pts));
  newpts(:,1) = r.*cos(dev) + pivot(1);   % New rectangular coordinates,
  newpts(:,2) = r.*sin(dev) + pivot(2);   %   restoring pivot
  i = find(~isfinite(rowsum(newpts)));
  if (~isempty(i))
    newpts(i,:) = ones(length(i),1)*pivot;
  end;
  
%   newpts = round(newpts*1e6)/1e6;

  if (doplot)
    scatter(savepts);
    axis equal;
    puttitle('Before rotation');
    
    scatter(newpts);
    axis equal;
    puttitle('After rotation');
  end;
  
  return;

