% DIGITIZE: Get point coordinates from a graphics image read from a graphics 
%           file.  Coordinates are scaled to a max range of 0-1.
%
%     Syntax: crds = digitize(filename,{fmt})
%
%           filename - name of graphics file, in quotes, including full path.
%           fmt =      optional file format, in quotes [default = 'tif'].
%                        Allowable: 'jpg','jpeg','tif','tiff','bmp','png',
%                        'hdf','pcx','xwd'.
%           -----------------------------------------------------------------
%           crds -     [n x 2] matrix of point coordinates.
%

% RE Strauss, 10/5/96
%   3/12/01 - Major rewrite using imread function.
%   4/1/02 -  Delete unneeded statements;
%             flip y-axis.


function crds = digitize(filename,fmt)
  if (nargin < 2) fmt = []; end;

  if (isempty(fmt))
    fmt = 'tif';
  end;

  map = imread(filename,fmt);
  colormap((1:256)' * [1 1 1]/256);
  image(map);
  axis equal;
  hold on;

  disp('  Digitize data points, press any key to stop');

  crds = [];
  while (1)                             % Capture points
    [x,y,button] = ginput(1);
    if (isempty(button))                % End of points
      break;
    end;
      crds = [crds; x y];                 % Accumulate point coordinates
      plot(x,y,'+r');
  end;

  crds(:,2) = -crds(:,2);               % Flip y-axis
  [N,P] = size(crds);
  smallest = min(crds);
  largest =  max(crds);
  range = max(largest-smallest);
  for pt = 1:N
    crds(pt,:) = (crds(pt,:)-smallest)./range;
  end;

  return;
