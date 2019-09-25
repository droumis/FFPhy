% REGROT: Given a configuration of points, of which two are landmarks,
%         registers the configuration so that the first landmark is located at 
%         (0,0), and rotates the configuration so that the second landmark is 
%         horizontal to the first.
%
%     Usage: [newpts,theta] = regrot(p1,p2,pts,{doplot})
%
%         p1 =     index to registration landmark within configuration.
%         p2 =     index to rotation landmark.
%         pts =    [N x 2] matrix of coordinates of point configuration.
%         doplot = optional boolean variable indicating, if true, that
%                    plots of the point configuration before and after
%                    rotation are to be produced [default = 0].
%         --------------------------------------------------------------
%         newpts = [N x 2] matrix of registered & rotated points.
%         theta =  angle of rotation.
%

% RE Strauss, 6/26/98
%   3/14/02 - return angle of rotation.
%   5/6/03 -  added error message and optional plots.

function [newpts,thetarot] = regrot(p1,p2,pts,doplot)
  if (~nargin) help regrot; return; end;
  
  if (nargin < 4) doplot = []; end;
  
  if (isempty(doplot)) doplot = 0; end;
  
  [N,P] = size(pts);
  if (P~=2)
    error('  RegRot: 2-dimensional point configurations only.');
  end;

  if (doplot)
    figure;
    plot(pts(:,1),pts(:,2),'ko',...
         pts(p1,1),pts(p1,2),'r*',...
         pts(p2,1),pts(p2,2),'k*');
    puttitle('Before rotation');
    putbnds(pts(:,1),pts(:,2));
    axis equal;
  end;

  pts = pts - ones(N,1)*pts(p1,:);            % Register on first landmark

  theta = anglerotation([0 0],[1 0],pts,1);   % Current angles with horizontal
  thetarot = -theta(p2);                      % Angle of rotation to second landmark
  theta = theta - theta(p2);                  % Rotate to second landmark

  r = sqrt(pts(:,1).^2 + pts(:,2).^2);        % Distances from origin

  newpts = pts;
  newpts(:,1) = r.*cos(theta);                % Rectangular coordinates
  newpts(:,2) = r.*sin(theta);
  newpts(p1,:) = [0 0];
  
  if (doplot)
    figure;
    plot(newpts(:,1),newpts(:,2),'ko',...
         newpts(p1,1),newpts(p1,2),'r*',...
         newpts(p2,1),newpts(p2,2),'k*');
    puttitle('After rotation');
    putbnds(newpts(:,1),newpts(:,2));
    axis equal;
  end;

  return;

