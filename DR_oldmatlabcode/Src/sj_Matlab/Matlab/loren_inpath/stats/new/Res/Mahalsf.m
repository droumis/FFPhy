% MAHALSF: Finds size-free Mahalanobis distances
%
%     Usage:  D2sf = mahalsf(X,grps,{kindsize})
%
%        X =    [n x p] data matrix (obs x vars).
%        grps = row or column vector of group identifiers.
%        kindsize =    kind of size vector: 'w' for within-group or
%                        or 'a' for among-group [default = 'w'].
%        ---------------------------------------------------------
%        D2sf = [k x k] matrix of size-free Mahalanobis distances.
%

% RE Strauss, 4/6/01
%   2/19/02 - added 'kindsize' option.

function D2sf = mahalsf(X,grps,kindsize)
  if (nargin < 3) kindsize = []; end;
  
  [loadings,percvar,scores,fscores,D2sf] = sizefree(X,grps,[],[],kindsize);

  return;

