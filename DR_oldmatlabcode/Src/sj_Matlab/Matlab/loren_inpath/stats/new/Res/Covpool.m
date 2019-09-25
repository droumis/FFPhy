% COVPOOL:  Calculates pooled within-group covariance matrix.  
%           Missing values are replaced by within-group means.
%
%     Usage: [Cp,mean_W,issing] = covpool(X,grps,{varonly})
%
%         X =       [n x p] data matrix.
%         grps =    [n x 1] vector of group identifiers for k groups.
%         varonly = optional boolean flag indicating, if true, that only pooled 
%                     within-group variances are to be returned on the diagonal
%                     of Cp [default = 0].
%         ---------------------------------------------------------------------
%         Cp =      [p x p] pooled within-group covariance matrix.
%         mean_W =  [k x p] matrix of within-group means.
%         issing =  boolean flag indicating, if true, that pooled covariance 
%                     matrix is singular.
%

% RE Strauss, 5/26/99
%   11/29/99 - changed calling sequence.
%   5/4/00 -   added option of variances-only.
%   12/12/00 - added check for singularity and return of within-group means.

function [Cp,mean_W,issing] = covpool(X,grps,varonly)
  if (nargin < 3) varonly = []; end;

  get_issing = 0;
  if (nargout > 2)
    get_issing = 1;
  end;

  if (isempty(varonly))
    varonly = 0;
  end;

  index = uniquef(grps);
  ngrps = length(index);            % Number of groups
  [nobs,nvars] = size(X);           % Numbers of observations & variables
  ndists = ngrps*(ngrps-1)/2;       % Pairwise combinations

  N =  zeros(ngrps,1);              % Within-group sample sizes
  Cp = zeros(nvars,nvars);          % Pooled covariance matrix

  G = design(grps);                 % Design matrix
  mean_W = (G'*G)\G'*X;             % Within-group means
  for g = 1:ngrps                   % Pooled within-group covariance matrix
    index = find(G(:,g));             % Indices to nonzero elements
    n = length(index);                % Sample size of current group
    N(g) = n;                         % Stash sample size
    Y = X(index,:);                   % Get data for current group

    [i,j] = find(~finite(Y));         % Handle missing data
    if (~isempty(i))
      mean_W(g,:) = means(Y);
      for k = 1:length(i)
        Y(i(k),j(k)) = mean_W(g,j(k));
      end;
    end;

    Ydev =  Y - ones(n,1)*mean(Y);    % Deviations
    Cp = Cp + Ydev'*Ydev;             % Augment within-group sum-of-squares
  end;
  Cp = Cp/(nobs-ngrps);             % Final pooled covariance matrix

  if (varonly)                      % Reduce only to variances on diagonal
    Cp = diag(diag(Cp));
  end;

  if (get_issing)                   % Check for singular matrix
    sizeCp = size(Cp,1);
    rankCp = rank(Cp);
    issing = 0;
    if (rank(Cp) < size(Cp,1))        % If covar matrix is singular,
      issing = 1;                     %   set flag
    end;
  end;

  return;

