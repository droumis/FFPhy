% PLOTGRPS: Produces a plot for several groups.
%
%     Syntax: plotgrps(X,Y,grps,{grpid},{symbol},{options})
%
%           X,Y  =    abscissa and ordinate vectors of length n.
%           grps =    grouping (classification) variable (vector of length n).
%                       Groups need not be sequenced or segregated.  Can be a 
%                       null vector of a single group.
%           grpid =   optional matrix of group identifier (1 row per grp) for 
%                       plotting, in ascii collating sequence.
%           symbol =  optional vector of symbols (either text strings or numeric 
%                       integers, 1 per grp) for plotting.  The symbols must appear
%                       in the ascii collating sequence of the grouping vector.  
%           options = optional vector of boolean flags specifying options,
%                       in any combination.  Presence of any 'options' vector
%                       resets the defaults.
%                         1) plot points [default]
%                         2) plot convex hulls [default]
%                         3) print group labels [default]
%                         4) plot centroids with '+' symbol [default]
%                         5) plot mediods with 'x' symbol
%                         6) print centroid or mediod labels
%                         7) plot major axes
%                         8) plot predictive regression lines
%                         9) plot 95% confidence ellipses of centroids
%                        10) plot 95% confidence ellipses of data
%                        11) plot 95% standard-error bars parallel to axes
%                        12) square plot
%                        13) equal axis units
%

% RE Strauss
%   5/31/99 - major general revision.
%  11/29/99 - changed calling sequence.
%   4/ 5/00 - corrected problem with hull for single point.
%  10/13/00 - plot groups in collating sequence;
%             check for agreement between group-label and group-id matrices.
%  10/16/00 - allow for row vectors for grp-id and symbol matrices.

function plotgrps(X,Y,grps,grpid,symbol,options)
  if (nargin < 4) grpid = []; end;
  if (nargin < 5) symbol = []; end;
  if (nargin < 6) options = []; end;

  default_options = [1 1 1 1];          % Default options
                                        % Set option flags
  [plot_points,plot_hulls,grp_labels,plot_centroids,plot_mediods,centr_labels, ...
    plot_majaxes,plot_regr,plot_cis,plot_ellipses,plot_se_bars,...
    square_plot,equal_axes] = pltgrpop(options,default_options);     

  if (isempty(grps))                    % If no groups, 
    grps = ones(size(X));               %   initialize group vector
  end;

  n = length(grps);
  if (length(X)~=n | length(Y)~=n)
    error('  PLOTGRPS: group and data vectors must be of same length.');
  end;

  indx = (finite(X) & finite(Y));       % Remove NaN's
  X = X(indx);
  Y = Y(indx);
  grps = grps(indx);

  if (isvector(grpid))                  % Convert row vector of grp id's to col vector
    grpid = grpid(:);
  end;
  if (isvector(symbol))
    symbol = symbol(:);
  end;
  if (isempty(grpid) & ~isempty(symbol))
    grpid = symbol;
  end;

  XY = [min(X) min(Y); max(X) max(Y)];  % Axis ranges

  [indx,freq] = uniquef(grps,1);        % Find groups and frequencies,
  ngrps = length(indx);                 %   in collating sequence

  if (~isempty(grpid))
    if (size(grpid,1) ~= ngrps)
      grpid
      ngrps
      error('  PLOTGRPS: size of group-label matrix must match number of groups.');
    end;
  end;

  hold on;
  for g = 1:ngrps                       % Cycle thru groups
    pts = find(grps == indx(g));          % Isolate points for current group
    x = X(pts);
    y = Y(pts);

    if (plot_points)
      pltgrppp(x,y,g,symbol);             % Plot data points
    end;

    if (plot_hulls)                       % Plot convex hulls
      deltax = 0.02 * range(X);
      pltgrpch(x,y,g,grp_labels,grpid,plot_points,symbol,deltax);
    end;

    if (plot_centroids)                   % Plot group centroids
      pltgrpcn(x,y,g,centr_labels,grpid,range([X Y]));
    end;

    if (plot_mediods)                     % Plot group mediods
      pltgrpmd(x,y,g,centr_labels,grpid,range([X Y]));
    end;

    if (freq(g) > 2)
      if (plot_regr)                      % Plot regression lines
        XY = pltgrppr(x,y,XY);
      end;
      
      if (plot_majaxes | plot_cis | plot_ellipses)  % Major axes or ellipses
        XY =pltgrpma(x,y,XY,plot_majaxes,plot_cis,plot_ellipses, ...
                      plot_regr,g,plot_hulls,grp_labels,grpid);
      end;

      if (plot_se_bars)                   % Plot standard-error bars
        pltgrpse(x,y,g,plot_centroids,centr_labels,grpid,range([X Y]));
      end;

    end;  % if more than 2 points
  end;  % Next group

  putbnd(XY);                           % Set plot axes
  set(gca,'box','on');

  if (square_plot)                      % Other plot options
    axis('square');
  end;
  if (equal_axes)
    axis('equal');
  end;

  hold off;                             % End plot

  return;
