% MAHALF: Objective function for finding Mahalanobis distances among groups.
%         Note: D2 is summed across variables, and therefore is dependent on
%         (not adjusted for) the number of variables.
%
%     Syntax: D2 = mahalf(X,grps,not_used1,not_used2)
%
%        grps =     row or column vector of group identifiers.
%        X =        [n x p] data matrix (obs x vars).
%        not_used = extra parameters passed by bootstrp().
%        ---------------------------------------------------------------------
%        D2 =       [1 x g*(g-1)/2] row vector representing the upper triangular
%                     off-diagonal portion of the [g x g] symmetric distance matrix.
%

% RE Strauss, 3/20/98
%   11/27/99 - handle singular pooled within-grp covar matrix.
%   11/16/00 - added warning message for singular covar matrix.
%   12/12/00 - use covpool() to find pooled within-group covar matrix.

function D2 = mahalf(X,grps,not_used1,not_used2)
  index = uniquef(grps);
  ngrps = length(index);                  % Number of groups

  D2 = zeros(1,ngrps*(ngrps-1)/2);        % Lower off-diagonals of distance matrix

  [Cp,mean_W,is_sing] = covpool(X,grps);  % Pooled within-group covar matrix
                                          %   and within-group means

  if (is_sing)                            % If covar matrix is singular,
    disp('  MAHALF warning: covariance matrix is singular;');
    disp('     reducing data matrix to best subset of variables.');

    v = steprank(Cp);                       % Get best subset of variables
    X = X(:,v);                             % Reduce data to subset

    [Cp,mean_W] = covpool(X,grps);          % Reduced pooled covar matrix
  end;

  i = 1;
  for g1 = 1:(ngrps-1)                    % Stash distances
    for g2 = (g1+1):ngrps
      diff = mean_W(g1,:) - mean_W(g2,:);
      D2(i) = diff * (Cp \ diff');
      i = i+1;
    end;
  end;

  return;

