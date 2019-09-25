% CLASSIFK: Classifies each of the "known" observations by jackknifing it and 
%           determining its group membership based on Mahalanobis distances.  
%           Called by classify().
%
%     Usage: results = classifk(X,grps,grpvals,sizeadj)
%
%           X =         [n x p] data matrix (obs x vars) for the n "known" training 
%                         observations.
%           grps =      [n x 1] column vector of group identifiers for k groups.
%           grpvals =   [k x 1] vector of unique group identifiers in 'grps'.
%           sizeadj =   boolean flag indicating, if true, that Mahalanobis 
%                         distances are to be based on residuals from a "size-free" 
%                         discriminant analysis [default = 0].
%           -----------------------------------------------------------------------
%           results =   [n x k] boolean matrix specifying the classification of 
%                         each observation into one of the k groups, allowing 
%                         for possible ties.
%

% RE Strauss, 7/13/98
%   11/29/99 - changed calling sequence.

function results = classifk(X,grps,grpvals,sizeadj)
  nobs = size(X,1);                       % Number of observations
  ngrps = length(grpvals);                % Number of groups
  newval = max(grps)+1;                   % Group identifier for "unknown"

  results = zeros(nobs,ngrps);            % Allocate output matrix

  for i = 1:nobs                          % Cycle thru observations
    x = [X; X(i,:)];                        % Move current obs to end
    x(i,:) = [];

    g = grps;                               % Assign to its own grp
    g(i) = [];
    g = [g; newval];

    if (sizeadj)                            % Get Mahalanobis distances
      [l,p,s,f,D2] = sizefree(x,g,1);
    else
      D2 = mahal(x,g);
    end;

    [minD2,closest] = min(D2(ngrps+1,1:ngrps));   % Set boolean classification
    results(i,closest) = ones(1,length(closest)); %   flag(s)
  end;

  return;

