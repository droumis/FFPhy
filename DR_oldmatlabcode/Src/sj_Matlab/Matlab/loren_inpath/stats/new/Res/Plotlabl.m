% Plotlabl: Modification of plot() command to plot the labels provided as rows  
%           of a corresponding matrix.  If only a single label value is provided, 
%           it is used for all points.  Default fontsize is 10.
%
%     Syntax:  plotlabl(x,y,labl,{fontsize})  OR
%              plotlabl([x,y],labl,{fontsize})
%

% RE Strauss, 10/10/97
%  9/7/99 -   miscellaneous changes for Matlab v5.
%  12/4/01 -  removed clf before plot statement.
%  9/24/02 -  remove plotting of four corners to initialize plot.
%  10/22/02 - change axis rescaling to use 'putbnds'.

function plotnum(x,y,labl,fontsize)
  if (nargin < 3) labl = []; end;
  if (nargin < 4) fontsize = []; end;
  
  if (isvector(x))
    x = x(:);
    if (isvector(y))
      y = y(:);
    end;
  end;

  [xr,c] = size(x);
  if (c>1)
    fontsize = labl;
    labl = y;
    y = x(:,2);
    x = x(:,1);
  end;

  if (isempty(fontsize))
    fontsize = 10;
  end;

  if (length(x)~=length(y))
    error('  PLOTLABL: coordinate vectors must be of same length');
  end;

  if (~ischar(labl))                      % Convert numeric to text labels
    labl = tostr(labl);
  end;

  [lr,c] = size(labl);
  if (lr~=1 & lr~=xr)
    error('  PLOTLABL: error in length of label matrix');
  end;

  % Expand single label to vector of identical labels
  if (lr==1)
    labl = char(ones(length(x),1)*double(labl));
  end;

  % Produce plot to get axis bounds
  plot(x,y,'w.');
  putbnds(x,y);
  v = axis;
  deltax = 0.011 * (v(2)-v(1));

  % Plot labels
  
  hold on;
  for i = 1:length(x)
    lab = labl(i,:);
    if (~ischar(lab))
      lab = int2str(round(lab));
    end;
    h = text(x(i)-deltax,y(i),lab);
    set(h,'fontsize',fontsize);
  end;
  hold off;

  return;
