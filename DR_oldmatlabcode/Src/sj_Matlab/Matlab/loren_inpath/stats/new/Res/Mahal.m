% MAHAL: Mahalanobis distances among groups, based on the pooled covariance  
%        matrix among all groups.  Reduces the number of variables to an optimal 
%        subset if the pooled covariance matrix is singular.
%        Note: D2 is summed across variables, and therefore is dependent on
%          (not adjusted for) the number of variables.
%
%     Syntax: [D2,CI,prob,power] = mahal(X,grps,{iter},{CI_level})
%
%        X =         [n x p] data matrix (obs x vars).
%        grps =      row or column vector of group identifiers.
%        iter =      number of bootstrap and randomization iterations [default=0].
%        CI_level =  percent width of confidence intervals [default=95].
%        ---------------------------------------------------------------------
%        D2 =    [g x g] symmetric distance matrix.
%        CI =    [g x g] matrix of low (lower triangular matrix) and
%                  high (upper triangular matrix) confidence limits.
%                  If iter=0, they are asymptotic estimates using a
%                    noncentral F-distribution.
%                  If iter>0, they are unbiased bootstrapped estimates.
%        prob =  [g x g] matrix of probabilities of H0:D2=0.
%                  If iter=0, they are asymptotic estimates based on an
%                    F-statistic.
%                  If iter>0, both asymptotic (lower triangular matrix) and
%                    randomized (upper triangular matrix) estimates are given.
%        power = [g x g] symmetric matrix of power estimates.
%

% Q: Because the pooled covar matrix is diminished as a function of the number
%    of groups, should it be estimated separately for each pair of groups?
% Tentative answer: No, because then the pairwise distances within the matrix 
%    will not be comparable.

% RE Strauss, 5/21/95
%   3/20/98 -   major rewrite.
%   11/27/99 -  handle singular covar matrix.
%   11/29/99 -  reversed X and grps in calling sequence.
%   6/13/00 -   added check for missing data.

function [D2,CI,prob,power] = mahal(X,grps,iter,CI_level)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  index = uniquef(grps);
  ngrps = length(index);                % Number of groups
  [nobs,nvars] = size(X);               % Numbers of observations & variables
  ndists = ngrps*(ngrps-1)/2;           % Pairwise combinations

  if (length(grps) ~= nobs)
    error('  MAHAL: Group vector and data matrix incompatible.');
  end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  else
    if (CI_level > 1)
      CI_level = 0.01 * CI_level;
    end;
  end;

  getCI =     0;                        % Initialize action flags
  getprob =   0;
  getpower =  0;
  bootstrap = 0;
  randomize = 0;

  if (nargout > 1)                      % Set action flags
    getCI = 1;
  end;
  if (nargout > 2)
    getprob = 1;
  end;
  if (iter > 0)
    if (getCI | getprob)
      bootstrap = 1;
    end;
    if (nargout > 3)
      getpower = 1;
    end;
  end;

  if (misscheck(X,grps))
    error('  MAHAL: data matrix or grouping vector contains missing data.');
  end;

  % ===== Mahalanobis distance matrix ===== %

  Dvect = mahalf(X,grps);           % Triangular matrix
  D2 = trisqmat(Dvect);             % Distance matrix

  % ===== Asymptotic probabilities and confidence intervals ===== %

  if (getprob & ~bootstrap)
    disp('  Asymptotic probabilities...');
    prob = mahalpr(X,grps,D2);
  end;

  if (getCI & ~bootstrap)
    disp('  Asymptotic confidence intervals...');
    CI = mahalci(X,grps,D2,CI_level);
  end;

  % ===== Bootstrap ===== %

  if (bootstrap)
    procs = [getCI, getprob, getpower, 0];
    alpha = 1-CI_level;
    nulldist = 0;
    [CI,prob,power]= bootstrp('mahalf',procs,iter,alpha,X,grps,nulldist);
  end;

  if (getCI)                 % Put results into square sym-matrix form
    CI = trisqmat(CI');
  end;
  if (getprob & bootstrap)
    prob = trisqmat(prob');
  end;
  if (getpower)
    power = trisqmat(power');
  end;

  return;
