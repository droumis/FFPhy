% PATHNN:   Nearest-neighbor path among a set of n points in p dimensions, 
%           given either the point coordinates (for which euclidean distances 
%           are calculated) or a distance matrix.  The start point can be 
%           specified, but defaults to the first point.
%
%     Usage: [path,pathlen] = pathnn(crds | dist,{sn},{getall})
%
%           crds =    [n x p] matrix of point coordinates.
%             or
%           dist =    [n x n] matrix of interpoint distances.
%           sn =      optional starting point number [default = 1].
%           getall =  optional boolean flag indicating, if true, that all 
%                       non-unique NN paths are to be returned [default = 0 = 
%                       return only first path found].
%           -------------------------------------------------------------------
%           path =    row vector of indices of points in connected sequence.
%           pathlen = total path length.
%

% RE Strauss, 9/27/00
%   10/3/00 - return non-unique paths for discrete points; 
%               allow input of distance matrix or point coordinates.
%   6/5/03 -  save only one path if all paths not requested.

function [path,pathlen] = pathnn(crds,sn,getall)
  if (nargin < 2) sn = []; end;
  if (nargin < 3) getall = []; end;

  if (isempty(sn))     sn = 1; end;
  if (isempty(getall)) getall = 0; end;

  [n,p] = size(crds);
  if (sn<1 | sn>n | ~isintegr(sn))
    error('  PATHNN: invalid starting point');
  end;

  if (issqsym(crds))                      % If input matrix is square-symmetric
    dist = crds;                            % Distance matrix has been passed
    tol = 1e-6 * max(trilow(dist));
  else                                    % Else coordinates have been passed
    dist = eucl(crds);                      % Get euclidean distances among pts
    tol = 1e-6 * min(range(crds));          % Tolerance for judging equal nn dists
  end;

  dist = putdiag(dist,2*max(max(dist)));  % Put large values on diagonal

  path = zeros(1,n);                      % Begin to build a single path
  path(1) = sn;                           % First point of path is starting pt
  pathlen = 0;                            % Corresponding path lengths
  npath = 1;                              % Number of paths being examined

  for ipt = 2:n                           % Cycle thru remaining points  
    ip = 0;
    while (ip<npath)                        % Cycle thru paths
      ip = ip+1;

      p = 1:n;                                % Disregard pts already in path
      k = path(ip,(path(ip,:)>0));
      p(k) = [];

      d = dist(path(ip,ipt-1),p);             % Dists from last pt to others
      dmin = min(d);                          % NN dist
      i = find(abs(d-dmin)<tol);              % Indices of pts having that NN dist
      path(ip,ipt) = p(i(1));                 % Add first NN point to path
      pathlen(ip) = pathlen(ip) + d(i(1));    % Update path length
      
      if (length(i)>1 & getall)               % If NN dist not unique,
        for j = 2:length(i)                     % For other pts with identical NN dists,
          path = [path(1:ip,:); path(ip,:); path(((ip+1):npath),:)];        % Duplicate path
          pathlen = [pathlen(1:ip); pathlen(ip); pathlen(((ip+1):npath))];  %   and path len

          ip = ip+1;
          path(ip,ipt) = p(i(j));                       % Replace pt identifier 
          pathlen(ip) = pathlen(ip-1)-d(i(1))+d(i(j));  % Adjust path length

          npath = npath+1;
        end;
      end;
    end;
  end;

  if (npath>1)                          % If path is not unique,
    plen = pathlen;                       % Find shortest path of set
    [pathlen,i] = min(pathlen);
    if (getall)                           % If requested to return all shortest paths,
      i = find(plen==pathlen);              % Find them
    end;
    path = path(i,:);                     % Return shortest path(s)
  end;

  return;
