% GETPOLYG: Selects a series of points from the graph window using the
%           mouse.  Points can be grouped into polygons, with consecutive
%           points joined by line segments, or alternately digitized as a
%           point string.  When the right mouse button is pressed, the
%           polygon is closed and a [0,0] delimiter is added to the list
%           of coordinates.  Press the <RETURN> (=<ENTER>) key to end the
%           accumulation of points.  Point coordinates are proportionately
%           rescaled to a [0,1] interval across groups.
%
%     Usage: [crds,grps] = getpolyg
%
%           crds = [n x 2] matrix of standardized point coordinates
%           grps = [n x 1] vector of matching group identifiers
%

% RE Strauss, 1/18/95
%   8/20/99 - changed plot colors and button/key use for consistency with Matlab v5.

function [crds,grps] = getpolyg
  disp(sprintf('  Press: right mouse button to close polygon,'));
  disp(sprintf('         <Enter> key (with cursor in window) to end input.'));

  crds = [];                          % Return matrices
  grps = [];
  dupl_crds = [];                     % Duplicate crds w/o delimiters
  delim = [0 0];                        % Delimiter

  RIGHT = 2;                            % Right mouse button
  TRUE  = 1;
  FALSE = 0;

  startpolygon = TRUE;
  grp_id = 1;

  clf;                                  % Put graph on screen
  axis([0 1 0 1]);
  plot([0 0 1 1 0],[0 1 1 0 0],'-k');
  axis('off');
  hold on;

  while (TRUE)                          % Capture points
    [x,y,button] = ginput(1);

    if (isempty(button))                % End of points
      break;
    end;

    if (startpolygon)                   % First point
      xprev = x;
      yprev = y;
      xpolyfirst = x;                   % Begin polygon
      ypolyfirst = y;
      startpolygon = FALSE;
    end;

    if (button > 1)                     % End of current polygon
      x = xpolyfirst;
      y = ypolyfirst;
      startpolygon = TRUE;
    end;

    crds = [crds; x y];             % Accumulate point coordinates
    dupl_crds = [dupl_crds; x y];   % Duplicate for standardization
    grps = [grps; grp_id];
    if (startpolygon)                   % Delimiter between polygons
      crds = [crds; delim];
      grps = [grps; 0];
      grp_id = grp_id + 1;
    end;

    plot(x,y,'+r');
    plot([xprev;x],[yprev;y],'-k');
    xprev = x;
    yprev = y;
  end;

  [N,P] = size(crds);
  if (startpolygon)                     % Delete final delimiter
    crds(N,:) = [];
    grps(N) = [];
    N = N-1;
  end;

  smallest = min(dupl_crds);
  largest =  max(dupl_crds);
  range = max(largest-smallest);
  for pt = 1:N
    if (~all(crds(pt,:) == delim))
      crds(pt,:) = (crds(pt,:)-smallest)/range;
    end;
  end;

  return;

