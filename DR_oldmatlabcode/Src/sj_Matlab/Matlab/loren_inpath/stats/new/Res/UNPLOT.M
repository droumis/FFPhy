% UNPLOT: Given axis coordinates and values and the coordinates of a set of
%         data points digitized from an image of a scatterplot, recovers the
%         data values.  Use [x,y]=impixel to get coordinates from the current image.
%
%     Usage: newcrds = unplot(axis_vals,axis_crds,crds)
%
%         axis_vals = vector of [xmin,xmax,ymin,ymax] reference values from axes.
%         axis_crds = [4 x 2] matrix of coordinates of xmin,xmax,ymin,ymax.
%         crds =      [n x 2] matrix of coordinates of data points.
%         -----------------------------------------------------------------------
%         newcrds =   [n x 2] matrix of recovered x,y data values from plot.
%

% RE Strauss, 7/13/95

function newcrds = unplot(axis_vals,axis_crds,crds)

  % Coordinates of x-axis
  px = axis_crds(1,1);
  py = axis_crds(1,2);
  qx = axis_crds(2,1);
  qy = axis_crds(2,2);
  x_axis = [(py-qy),(qx-px),(px*qy-qx*py)];

  % Coordinates of y-axis
  px = axis_crds(3,1);
  py = axis_crds(3,2);
  qx = axis_crds(4,1);
  qy = axis_crds(4,2);
  y_axis = [(py-qy),(qx-px),(px*qy-qx*py)];

  [n,p] = size(crds);
  if (p~=2)
    error('  Coordinates matrix of wrong size or shape');
  end;

  x_proj = zeros(n,2);                % Projections onto x-axis
  y_proj = zeros(n,2);                % Projections onto y-axis

  for i = 1:n                         % Cycle thru data points
    px = crds(i,1);
    py = crds(i,2);

    % Line through data point orthogonal to x_axis
    x_orth = [-x_axis(2),x_axis(1),(x_axis(2)*px-x_axis(1)*py)];
    % Intersection of x_axis and orthogonal line, =x projection
    dlm = x_axis(1)*x_orth(2) - x_axis(2)*x_orth(1);
    x_proj(i,:) = [(x_axis(2)*x_orth(3)-x_axis(3)*x_orth(2))/dlm, ...
                   (x_axis(3)*x_orth(1)-x_axis(1)*x_orth(3))/dlm];

    % Line through data point orthogonal to y_axis
    y_orth = [-y_axis(2),y_axis(1),(y_axis(2)*px-y_axis(1)*py)];
    % Intersection of y_axis and orthogonal line, =y projection
    dlm = y_axis(1)*y_orth(2) - y_axis(2)*y_orth(1);
    y_proj(i,:) = [(y_axis(2)*y_orth(3)-y_axis(3)*y_orth(2))/dlm, ...
                   (y_axis(3)*y_orth(1)-y_axis(1)*y_orth(3))/dlm];
  end;

  x_length = eucl([axis_crds(1,1),axis_crds(1,2)], ...
                  [axis_crds(2,1),axis_crds(2,2)]);
  x_min = axis_vals(1);
  x_max = axis_vals(2);
  x_units = x_max - x_min;
  x = x_min + eucl(x_proj,[axis_crds(1,1),axis_crds(1,2)])*x_units/x_length;

  y_length = eucl([axis_crds(3,1),axis_crds(3,2)], ...
                  [axis_crds(4,1),axis_crds(4,2)]);
  y_min = axis_vals(3);
  y_max = axis_vals(4);
  y_units = y_max - y_min;
  y = y_min + eucl(y_proj,[axis_crds(3,1),axis_crds(3,2)])*y_units/y_length;

  newcrds = [x,y];

  return;
