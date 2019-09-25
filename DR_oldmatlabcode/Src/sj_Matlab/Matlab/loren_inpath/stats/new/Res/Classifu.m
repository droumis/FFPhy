% CLASSIFU: Classifies "unknown" observations, singly or in groups, based on a  
%           training set of observations.  Called by classify().
%
%     Usage: results = classifu(X,grps,Xu,ugrps,ntgrps,ugrpvals,sizeadj)
%           X =         [n x p] data matrix (obs x vars) for the n "known"  
%                         training observations.
%           grps =      [n x 1] vector of group identifiers for k training grps.
%           Xu =        [m x p] data matrix for the m "unknown" observations.
%           ugrps =     [m x 1] vector of group identifiers for "unknown"
%                         observations; if not passed, the unknowns are 
%                         classified individually.
%           ntgrps =    number of training groups.
%           ugrpvals =  vector of unique group identifiers in 'grps'.
%           sizeadj =   boolean flag indicating, if true, that Mahalanobis 
%                         distances are to be based on residuals from a 
%                         "size-free" discriminant analysis.
%           ---------------------------------------------------------------------
%           results =   [m x k] boolean matrix specifying the classification of 
%                         each "unknown" observation into one of the k training 
%                         groups, allowing for possible ties.
%

% RE Strauss, 7/13/98
%   11/29/99 - changed calling sequence.

function results = classifu(X,grps,Xu,ugrps,ntgrps,ugrpvals,sizeadj)
  nugrps = length(ugrpvals);              % Number of groups of unknowns
  nuobs = size(Xu,1);                     % Number of unknown observations
  results = zeros(nuobs,ntgrps);          % Allocate output matrix

  for iu = 1:nugrps                       % Cycle thru unknowns groups
    i = find(ugrps == ugrpvals(iu));        % Isolate unknowns for current grp

    g = [grps; ugrps(i)];                   % Append unknown(s) to training matrices
    x = [X; Xu(i,:)];

    if (sizeadj)                            % Get Mahalanobis distances
      [l,p,s,f,D2] = sizefree(x,g,1);
    else
      D2 = mahal(x,g);
    end;

    [minD2,closest] = min(D2(ntgrps+1,1:ntgrps));
    results(i,closest) = ones(length(i),length(closest));
  end;

  return;

