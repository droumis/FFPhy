% IMAGEBOX: Given an RGB image and the coordinates of the lower-left and upper-right
%           corners of a rectangle, returns (and optionally plots) the sub-image
%           corresponding to the rectangle.
%
%     Usage: boxrgb = imagebox(rgb_image,corners,{noflip},{doplot})
%
%           rgb_image = 3D RGB image matrix.
%           corners =   [2 x 2] matrix of the LL and UR rectangle corners.
%           noflip =    optional boolean flag indicating, if true, that the RGB image
%                         is not to be flipped upside down before applying coordinates
%                         [default = 0].
%           doplot =    optional boolean flag indicating that plot of submatrix is to
%                         be produced [default = 0].
%           --------------------------------------------------------------------------
%           boxrgb =    3D RGB submatrix.
%

% RE Strauss, 6/26/02

function boxrgb = imagebox(rgb_image,corners,noflip,doplot)
  if (nargin < 3) noflip = []; end;
  if (nargin < 4) doplot = []; end;
  
  if (isempty(noflip))
    noflip = 0;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;
  
  if (~noflip)                              % Flip image upside-down
    rgb_image = rgb_image(size(rgb_image,1):-1:1,:,:);
  end;
  
  xmin = min(corners(:,1));
  xmax = max(corners(:,1));
  ymin = min(corners(:,2));
  ymax = max(corners(:,2));
  
  boxrgb = rgb_image(ymin:ymax,xmin:xmax,:);
  
  if (doplot)
    figure;
    image(boxrgb);
  end;

  return;
  