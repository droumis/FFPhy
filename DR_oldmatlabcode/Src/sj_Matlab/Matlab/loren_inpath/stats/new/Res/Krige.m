% KRIGE: Use universal kriging to predict missing values from a spatial matrix.
%        Assumes linear drift within a minimally localized region of the map, and
%        a linear semivariogram.
%
%     Syntax: Z = krige(X,Y,Z,{nest},{ar})
%
%           X =    vector of x coordinates.
%           Y =    vector of y coordinates.
%           Z =    vector of corresponding z coordinates,
%                    with missing values to be kriged coded as NaN or inf.
%           nest = optional number of nearest points used to estimate missing 
%                    value [default = 5, the minimum].
%           ar =   optional vector describing the X,Y "aspect ratio" for scaling 
%                    the axes of the input matrices; axes will be scaled from 
%                    zero to these values; use [1 1] for equal scalings
%                    [default = no scaling].
%           --------------------------------------------------------------------
%           Z =    vector of z coordinates with missing values
%                    replaced by predicted values.
%

% Davis,JC. 1986. Statistics and Data Analysis in Geology, 2nd ed,
%   pp. 239-248, 383-405.  Wiley.

% RE Strauss, 10/30/95
%   9/20/99 - update handling of null input arguments.

function Z = krige(X,Y,Z,nest,ar)
  if (nargin < 4) nest = []; end;
  if (nargin < 5) ar = []; end;

  if (isempty(nest))
    nest = 5;                         % Nmbr nearest pts used to est missing value
  end;

  [r,c] = size(X);                    % Convert vectors to columns if necessary
  if (r==1 & c>1)
    X = X';
  end;
  [r,c] = size(Y);
  if (r==1 & c>1)
    Y = Y';
  end;
  [r,c] = size(Y);
  if (r==1 & c>1)
    Z = Z';
  end;

  if (length(ar)>0)                   % Scale X,Y axes if requested
    minX = min(X);
    minY = min(Y);
    X = ar(1)*(X-minX)/(max(X)-minX);
    Y = ar(2)*(Y-minY)/(max(Y)-minY);
  end;

  indx = find(~finite(Z));            % Locate missing Z values
  if (isempty(indx))
    error('KRIGE: No missing values to be kriged.');
  end;

  xobs = [X Y];                       % Matrix of observed (control) points
  xmiss = xobs(indx,:);               % Coordinates of points to be kriged
  xobs(indx,:) = [];                  % Remove missing points
  z = Z;
  z(indx) = [];
  n = size(xobs,1);                   % Number of observed points
  m = size(xmiss,1);                  % Number of missing points

  if (n < nest)
    error('KRIGE: Number of observed points < requested number of nearest points');
  end;

  eucldist = eucl(xmiss,xobs);        % Distances among missing and obs pts

  semivar = zeros(nest*(nest-1)/2,2); % Alloc semivariogram data matrix
  Gobs = zeros(nest,nest);            % Alloc gamma (semivariance) matrix for obs pts
  Gmiss = zeros(nest,1);              % Alloc gamma matrix for kriged pts
  singular = 0;

  for misspt = 1:m                    % Cycle thru missing points
    [distnear,id] = sort(eucldist(misspt,:)');  % Find set of nearest obs points
    nearest = id(1:nest);             % Indices of nearest points
    dist = eucl(xobs(nearest,:));     % Distances among nearest points

    k = 0;
    for i = 1:(nest-1)                % Pairwise semivariances among nearest pts
      for j = (i+1):nest
        k = k+1;
        semivar(k,1) = dist(i,j);
        semivar(k,2) = (z(nearest(i))-z(nearest(j))).^2/2;
      end;
    end;                              % Regress thru origin for linear semivariogram
    slope = (semivar(:,1)'*semivar(:,2)) / (semivar(:,1)'*semivar(:,1));

    for i = 1:(nest-1)                % Predicted semivariances for obs pts
      for j = (i+1):nest
        gamma = slope * dist(i,j);
        Gobs(i,j) = gamma;
        Gobs(j,i) = gamma;
      end;
    end;
    for i = 1:nest
      gamma = slope * eucldist(misspt,nearest(i));
      Gmiss(i) = gamma;
    end;
    xnear = xobs(nearest,:) - ones(nest,1)*xmiss(misspt,:);  % Centered coords

    A = [Gobs ones(nest,1) xnear; [ones(1,nest); xnear'] zeros(3,3)];  % A matrix
    B = [Gmiss; [1 0 0]'];                                             % B matrix

    if (rcond(A) > eps)               % If A matrix is conditioned (non-singular),
      W = A\B;                          % Solve for weights by Gaussian elimination
      w = W(1:nest);                    % Krige weights
    else                              % Else use weights proportional to dist^2
      distnear = 1./(distnear(1:nest).^2); % Col vector of recip squared dists
      w = distnear / sum(distnear);     % Reciprocal standardized
      singular = singular + 1;          % Singularity tally
    end;

    Z(indx(misspt)) = w'*z(nearest);  % Weighted predicted value for missing pt
  end;

  if (singular > 0)
    disp(sprintf('  %1.0f of %1.0f krigings were singular', singular,m));
  end;

  return;

