% CORRMAP:  Produces a color-patch representation of a correlation matrix,
%           where the fill density of each matrix cell is proportional to
%           the absolute value of the correlation coefficient.
%
%     Usage: corrmap(r)
%
%           r = square correlation matrix
%

% RE Strauss, 3/9/95

function corrmap(r)
  figure;                           % New figure
  n = max(size(r));                 % Size of correlation matrix
  colormap(bone);                   % Use density colormap
  map = flipud(colormap);           % Make black=max density rather than white
  colormap(map);
  m = size(map,1);                  % Rows in colormap
  c = ceil(m.*abs(r));              % Index into colormap matrix
  image(1:n,1:n,c);                 % Square plot%  axis('square');

  figure;                           % New figure
  r = flipud(r);                    % Flip matrix upside-down
  v = [.1:.2:.9];                   % Coutour every 0.2 units
  c = contour(r,v);                 % Contour plot

  return;

