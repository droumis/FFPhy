% GETAXES: Given an image specified by [X,map], allows the digitizing of the 
%          axis tick marks corresponding to [xlow,xhigh,ylow,yhigh].
%
%       Usage: ticks = getaxes(X,map)
%
%           X =     image matrix
%           map =   colormap for image
%           ticks = vector of axis tick marks corresponding to [xlow,xhigh,ylow,yhigh]
%

% RE Strauss, 8/12/95

function ticks = getaxes(X,map)
  sintheta = sin(pi/4);               % For rotations of 45 degrees
  costheta = cos(pi/4);

  disp('  Digitize xlow,xhigh, ylow,yhigh...');
  imshow(X,map);
  [x,y] = impixel;

  % Rotate entire configuration by 45 degrees

%  xin = x - min(x);
%  yin = -y - min(-y);
  xin = x;
  yin = -y;
  x = xin*costheta - yin*sintheta;
  y = yin*costheta + xin*sintheta;

  x = x(1:4);                         % Reduce to points of interest
  y = y(1:4);
  ticks = [x y];

  return;
