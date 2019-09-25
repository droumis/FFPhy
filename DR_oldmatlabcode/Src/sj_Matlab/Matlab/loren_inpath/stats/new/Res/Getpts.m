% GETPTS: Selects a series of points from the graph window using the
%         mouse.  Press the right mouse button to increment the group
%         identifier.  Press the <RETURN> (=<ENTER>) key to end the
%         accumulation of points.  Point coordinates are proportionately
%         rescaled to a [0,1] interval across groups.
%
%     Usage: [crds,grps] = getpts
%
%           crds = [n x 2] matrix of standardized point coordinates.
%           grps = [n x 1] corresponding group-identification vector.
%

% RE Strauss, 1/5/99, modified from getpolyg()
%   1/9/99 - added group-identifier output.
%   9/7/99 - miscellaneous changes for Matlab v5.

function [crds,grps] = getpts
  disp(sprintf('  Press <Enter> key (with cursor in window) to end input.'));

  crds = [];                            % Return matrices
  grps = [];
  dupl_crds = [];                       % Duplicate crds w/o delimiters

  clf;                                  % Put graph on screen
  axis([0 1 0 1]);
  plot([0 0 1 1 0],[0 1 1 0 0],'-w');
  axis('off');
  hold on;

  grp_id = 1;
  while (1)                             % Capture points
    [x,y,button] = ginput(1);
    if (isempty(button))                % End of points
      break;
    end;
    if (button > 1)                     % Increment group identifier
      grp_id = grp_id+1;
    else
      crds = [crds; x y];                 % Accumulate point coordinates
      grps = [grps; grp_id];
      dupl_crds = [dupl_crds; x y];       % Duplicate for standardization
      plot(x,y,'+r');
    end;
  end;

  [N,P] = size(crds);
  smallest = min(dupl_crds);
  largest =  max(dupl_crds);
  range = max(largest-smallest);
  for pt = 1:N
    crds(pt,:) = (crds(pt,:)-smallest)./range;
  end;

  return;

