% MOMENTS:  Given a matrix X, returns a matrix of column vectors containing
%           the r=1...R moments about the mean (or other specified point), for 
%           each column of X.  Moments for r>2 may optionally be scaled by the 
%           second moment (i.e., be made dimensionless).
%
%     Usage: mom = moments(X,{R},{dimless},{pts})
%
%           X =       [n x p] matrix.
%           R =       optional max number of moments [default=4].
%           dimless = optional boolean flag indicating, if true, that moments 
%                       3..R are to be scaled by s^r [default=0].
%           pts =     optional vector (length p) of points about which moments 
%                       are to be calculated [default = column means].  If a 
%                       constant is specified, it is used for all columns.
%           ------------------------------------------------------------------
%           mom =     [R x p] matrix of moments.
%

% RE Strauss, 4/20/95
%   11/23/99 - allow dimensionless moments.

function mom = moments(X,R,dimless,pts)
  if (nargin < 2) R = []; end;
  if (nargin < 3) dimless = []; end;
  if (nargin < 4) pts = []; end;

  if (isvector(X))                    % Transpose to col if row vector
    X = X(:);
  end;
  [nr,nc] = size(X);                  % Size of data matrix

  if (isempty(R))
    R = 4;
  end;
  if (isempty(dimless))
    dimless = 0;
  end;
  if (isempty(pts))
    pts = mean(X);
  end;

  if (~isvector(pts))
    if (isscalar(pts))
      pts = pts*ones(1,nc);
    else
      error('  MOMENTS: invalid vector of center points');
    end;
  end;

  dev = X - ones(nr,1)*pts;

  mom = zeros(R,nc);
  for r = 1:R
    mom(r,:) = sum(dev.^r)./nr;
    if (r>2 & dimless)
      mom(r,:) = mom(r,:)./(sqrt(mom(2,:)).^r);
    end;
  end;

  return;
