% PLOTDIST: Given a set of point coordinates and a distance-specification matrix, 
%           produces a plot of the points and interpoint distances.  Optionally 
%           plots specified values at the midpoints of the lines.
%
%     Usage: [xv,yv,v] = plotdist(crds,dists,{individ},...
%                                 {regpts},{title},{values},{digits},{fontsize})
%
%           crds =      [n x 2] matrix of point coordinates for n points for one 
%                         or more individuals.
%           dists =     [m x maxp] matrix of point-index specifications for m 
%                         pairwise interpoint distances.  If more than two points 
%                         are specified per row, the intermediate points are 
%                         assumed to be helping points.  
%                         Rows with <maxp points are buffered with zeros.
%           individ =   optional vector of indentifiers for individuals, 
%                         corresponding to 'crds'.  Separate plots are produced 
%                         for each individual, pausing after each one.
%           regpts =    optional vector of two point indices for registration 
%                         and rotation.  The point configuration is registered 
%                         on the first point and rotated so that the second point 
%                         is horizontal to the first [default = no action].
%           title =     optional character string for plot title.
%           values =    optional [m x 1] vector of numeric values for each 
%                         distance measure, or [n x c] character-string matrix.
%           digits =    number of significant digits for values [default = 3].
%           fontsize =  optional font size for printed values [default = 10].
%           ---------------------------------------------------------------------
%           xv,yv =     coordinates of points at which labels printed.
%           v =         values of labels printed.
%

% RE Strauss, 9/2/99
%   9/ 9/98 -   allow for helping points.
%   9/10/98 -   use linelabl() to print labels.
%   11/25/99 -  changed plot colors for Matlab v5.
%   6/1/01 -    add registration/rotation option.
%   8/1/01 -    added optional vector for individuals.

function [xv,yv,v] = plotdist(crds,dists,individ,regpts,title,values,digits,fontsize)
  if (nargin < 3) individ = []; end;
  if (nargin < 4) regpts = []; end;
  if (nargin < 5) title = []; end;
  if (nargin < 6) values = []; end;
  if (nargin < 7) digits = []; end;
  if (nargin < 8) fontsize = []; end;

  if (isempty(digits))                    % Default input argument values
    digits = 3;
  end;
  if (isempty(fontsize))
    fontsize = 10;
  end;

  plot_val = 0;
  v = [];

  if (~isempty(values))
    plot_val = 1;
    if (~ischar(values))                   % Convert labels to strings
      v = values;
      values = tostr(values,digits);
    end;
  end;
  
  if (~isempty(regpts))
    if (length(regpts(:))~=2)
      error('  PLOTDIST: Invalid registration-rotation point indices.');
    end;
  end;

  if (isempty(individ))
    individ = ones(size(crds,1));
  end;

  if (size(dists,2)==3)                   % If distance-specification matrix has
    dists = dists(:,2:3);                 %   3 cols, delete the first
  end;
  len_dists = size(dists,1);
  udist = uniquef(dists(:),1);            % Unique dists used in plot

  landmarks = [];
  help_pts = [];
  midpoints = zeros(len_dists,2);

  uindivid = uniquef(individ);
  
  for ui = 1:length(uindivid)             % Cycle thru individuals
    ic = find(individ==uindivid(ui));       % Isolate crds for individuals
    ucrds = crds(ic,:);
    if (~isempty(regpts))                   % Register and rotate
      ucrds = regrot(regpts(1),regpts(2),ucrds);
    end;

    for i = 1:len_dists                     % Tally lists of landmarks and helping points
      p = dists(i,:);
      pos = find(p>0);
      lenpos = length(pos);

      landmarks = [landmarks; p([pos(1) pos(lenpos)])'];
      if (lenpos > 2)
        help_pts = [help_pts; p(pos(2:lenpos-1))'];
        c = pathpts(ucrds(p(pos),:),3);
        midpoints(i,:) = c(2,:);
      else
        midpoints(i,:) = mean(ucrds(p(pos),:));
      end;
    end;
    landmarks = uniquef(landmarks,1);
    help_pts = uniquef(help_pts,1);

    xv = midpoints(:,1);
    yv = midpoints(:,2);

    x = ucrds(:,1);
    y = ucrds(:,2);

    figure;
    plotnum(x,y,fontsize);                  % Produce plot with distance identifiers
    putbnd(x,y);
    axis('equal');
    axis('off');

    deltax = 0.010 * (max(x)-min(x));       % Kerns for printed values
    deltay = 0.070 * (max(y)-min(y));

    figure;
    plot(ucrds(landmarks,1),ucrds(landmarks,2),'ko');
    hold on;
    plot(ucrds(help_pts,1),ucrds(help_pts,2),'kx');
    hold on;
    for i = 1:len_dists
      p = dists(i,:);
      pos = find(p>0);
      lenpos = length(pos);

      for j = 1:lenpos-1
        plot(ucrds(p([j,j+1]),1),ucrds(p([j,j+1]),2),'k');
      end;

      if (plot_val)
        endpoints = [ucrds(p(1),:); ucrds(p(lenpos),:)];
        label = values(i,:);
        position = [xv(i) yv(i)];

        linelabl(endpoints,label,position,deltax,deltay,fontsize);
      end;

    end;
    hold off;
    putbnd(x(udist),y(udist));
    axis('equal');
    axis('off');

    if (~isempty(title))
      puttitle(title);
    end;

    if (ui < length(uindivid))
      pause;
    end;
  end;

  return;

