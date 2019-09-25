% Putbnd3: Changes the [min,max] settings for the axes of a 3D plot to allow a 
%          buffer beyond the range of the data.
%
%     Usage: v = putbnd3(crds,{buffer},{nocall})
%
%         crds =   [n x 3] matrix of point coordinates used to create the plot.
%         buffer = optional buffer size, as proportion of ranges [default = 0.05].
%         nocall = optional boolean flag indicating, if true, that the axis 
%                    settings are to be returned but that the current plot is 
%                    to be left unaltered [default = 0].
%         ------------------------------------------------------------------------
%         v =      new axis ranges: [xmin xmax ymin ymax zmin zmax].
%

% RE Strauss, 6/6/03

function v = putbnd3(crds,buffer,nocall)
  if (nargin < 2) buffer = []; end;
  if (nargin < 3) nocall = []; end;
  
  if (isempty(buffer)) buffer = 0.05; end;
  if (isempty(nocall)) nocall = 0; end;
  
  x = crds(:,1);
  y = crds(:,2);
  z = crds(:,3);

  xmin = min(x);
  xmax = max(x);
  ymin = min(y);
  ymax = max(y);
  zmin = min(z);
  zmax = max(z);

  deltax = 0.05*(xmax-xmin);
  deltay = 0.05*(ymax-ymin);
  deltaz = 0.05*(zmax-zmin);

  xmin = xmin - deltax;
  xmax = xmax + deltax;
  ymin = ymin - deltay;
  ymax = ymax + deltay;
  zmin = zmin - deltaz;
  zmax = zmax + deltaz;
  
  v = [xmin xmax ymin ymax zmin zmax];
  if (~nocall)
    axis(v);
  end;

  return;
  