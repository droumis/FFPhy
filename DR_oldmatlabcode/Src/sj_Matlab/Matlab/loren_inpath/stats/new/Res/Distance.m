% DISTANCE: Converts a set of point coordinates to a set of Euclidean distances 
%           or square roots of areas, given the conversion specifications.
%
%       Usage: dist = distance(specs,crds,{individ},...
%                               {scalebar},{scaleunits},{missrefpt})
%
%           specs =     [m x 2], [m x 3] or [m x 4] set of interpoint distance 
%                         specifications:
%                           col 1: index of distance variable (1,2,...).
%                               2: index of first point.
%                               3: index of second point.
%                               4: indicator of type of output:
%                                   1 = cumulative interpoint distance.
%                                   2 = sqrt(area) of enclosed polygon.
%
%                         If col 1 is missing (i.e., if only 2 cols are present), 
%                           it is assumed to be a column vector 1:m.
%                         If col 4 is missing, it is assumed to be a column vector 
%                           of 1's (specifying all linear interpoint distances).
%           crds =      [n x p] set of point coordinates.
%           individ =   optional vector [length n] specifying k individuals, each 
%                         represented by the same number (np) of point coordinates 
%                         (so that np x k = n).  [Default: all coordinates for a  
%                         single individual].
%           scalebar =  optional vector of two point indices specifying 
%                         digitized scale bar.
%           scalelunits = optional units for scale bar [default = 1].
%           missrefpt = optional index specifying the point below which 
%                         'missing-data' points have been digitized.
%           ----------------------------------------------------------------------
%           dist =      [k x nd] row vector of distances.
%
%
%       Sample 'specs' matrix:    [1 1 3 1
%                                  2 4 5 1
%                                  2 5 6 1
%                                  3 3 4 2
%                                  3 4 5 2
%                                  3 5 8 2
%                                  3 8 6 2]
%       specifies 3 output variables:
%         1: the distance between points 1 and 3.
%         2: the distance between points 4 and 6, via point 5.
%         3: the sqrt(area) of the polygon enclosed by points 3, 4, 5, 8, and 6.
%

% RE Strauss, 4/5/97
%   9/19/99 -  updated handling of default input arguments.
%   12/1/99 -  allow for missing-data reference point and scale bar.
%   12/11/01 - fixed problem with missing points when have no scale bar.

function dist = distance(specs,crds,individ,scalebar,scaleunits,missrefpt)
  if (nargin < 3) individ = []; end;
  if (nargin < 4) scalebar = []; end;
  if (nargin < 5) scaleunits = []; end;
  if (nargin < 6) missrefpt = []; end;

  [m,c] = size(specs);
  [n,p] = size(crds);

  if (isempty(individ))                   % Coords for single individual
    individ = ones(n,1);
  else
    if (length(individ) ~= n)
      error('  DISTANCE: "crds" and "individ" matrices not compatible');
    end;
  end;

  rescale = 0;
  check_missing = 0;
  if (~isempty(scalebar))
    rescale = 1;
    if (max(scalebar) > n)
      error('  DISTANCE: scale bar point specification out of range');
    end;
    if (length(scalebar) ~= 2)
      error('  DISTANCE: scale bar must be vector of two point indices');
    end;
  end;
  if (~isempty(missrefpt))
    check_missing = 1;
    if (max(missrefpt) > n)
      error('  DISTANCE: scale bar point specification out of range');
    end;
  end;

  if (~isempty(scalebar) & isempty(scaleunits))
    scaleunits = 1;
  end;

  if (c~=4)
    if (c==2)
      specs = [[1:m]' specs];
      c = c+1;
    end;
    if (c==3)
      specs = [specs ones(m,1)];
    else
      error( '  DISTANCE: specification matrix must have 3 or 4 columns');
    end;
  end;

  distindx = specs(:,1);             % Extract info from specifications
  pts =  specs(:,2:3);
  disttype = specs(:,4);

  if (max(pts(:))>n | min(pts(:))<1) % Check point-spec range
    error('  DISTANCE: point specification out of range');
  end;

  if (max(disttype)<1 | min(disttype>2))
    error('  DISTANCE: invalid type specification');
  end;

  ndist = max(distindx);            % Number of output distance variables
  
  indlist = uniquef(individ);        % List of individual id's
  nind = length(indlist);           % Number of individuals

  dist = NaN*ones(nind,ndist);      % Allocate output distance matrix

  for id = 1:nind                   % Cycle thru individuals
    ind = indlist(id);
    crdind = crds(individ==ind,:);    % Extract points for current individual

    if (check_missing)                % Check for missing-data points
      i = find(crdind(:,2) < crdind(missrefpt,2));
      if (~isempty(i))                        % Replace pts below missing-data 
        if (rescale)
          j = find(i==scalebar(1) | i==scalebar(2));  % But not the scale bar pts
          if (~isempty(j))
            i(j) = [];
          end;
        end;
        if (~isempty(i))
          crdind(i,:) = NaN*ones(length(i),p);  %   reference pt with NaN's
        end;
      end;
    end;

    for v = 1:ndist                   % Cycle thru variables
      indx = find(distindx==v);         % Extract specs for current variable
      if (length(indx)>0)
        type = disttype(indx(1));         % Type of variable
        if (type==1)
          p1 = crdind(pts(indx,1),:);       % 1st endpoint
          p2 = crdind(pts(indx,2),:);       % 2nd endpoint
          delta = (p1-p2)';                 % Diffs along coordinate axes
          dist(id,v) = sum(sqrt(sum(delta.*delta)')'); % Euclidean distance
        elseif (type==2)
          p = pts(indx,:)';
          p = uniquef(p(:));
          dist(id,v) = sqrt(polyarea(crdind(p,:)));
        end;
      end;
    end;

    if (rescale)                      % Rescale with respect to scale bar
      scale_factor = eucl(crdind(scalebar,:)) ./ scaleunits;
      dist(id,:) = dist(id,:)./scale_factor;
    end;
  end;
  
  return;
