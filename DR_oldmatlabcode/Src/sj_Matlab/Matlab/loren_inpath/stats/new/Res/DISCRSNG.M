% DISCRSNG: Discriminant analysis for singular or near-singular conditions 
%           due to too many variables relative to the number of observations.
%           During each iteration the specified number of variables is randomly 
%           selected from the total available.
%     
%     Syntax:  [load,percvar,scores,CI_load,CI_percvar] ...
%                 = discrsng(X,grps,nvars,{ndf},{iter},{CI_level})
%
%        X =           [n x p] data matrix (obs x vars).
%        grps =        row or column vector of group identifiers.
%        nvars =       number of variables to be randomly selected from the 
%                        total available during each randomized replication.
%        ndf =         optional number of leading discriminant functions for
%                        which scores are desired (default = groups-1).
%        iter =        number of randomized iterations (default = 1).
%        CI_level =    percent width of confidence intervals (default = 95).
%        ---------------------------------------------------------------------
%        load =       [p x ndf] matrix of discriminant-function
%                        coefficients (columns) as vector correlations.
%        percvar =    [ndf x 1] vector of percents of total variance explained
%                        for discriminant functions.
%        scores =     [p x ndf] matrix of mean discriminant-function scores
%                        (columns), across randomized replications.
%        CI_load =    [p x 2*ndf] matrix of CI% confidence limits
%                        (asymmetric) of vector-correlation loadings, two columns
%                        per discriminant function (low, high).
%        CI_percvar = [p x 2] matrix of CI% confidence limits (asymmetric)
%                        of percents variance-explained.
%

% RE Strauss, 11/18/96
%   9/20/99 - update handling of null input arguments.
%  11/29/99 - changed calling sequence.

function [load,percvar,scores,CI_load,CI_percvar] = ...
            discrsng(X,grps,nvars,ndf,iter,CI_level)

  if (nargin < 3) ndf = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) CI_level = []; end;

  default_ndf = ngrps-1;          % Defaults for input arguments
  default_iter = 1000;
  default_CI_level  = 95;

  index = uniquef(grps);
  ngrps = length(index);          % Number of groups
  [nobs,nvars] = size(X);         % Numbers of observations & variables
  G = design(grps);               % ANOVA-type design matrix

  if (isempty(ndf))               % If ndf not passed,
    ndf = default_ndf;            %   set to default
  end;
  if (isempty(iter))              % If iter not passed, set to default
    iter = default_iter;
  end;
  if (isempty(CI_level))          % Check confidence level
    CI_level = default_CI_level;
  end;

  do_percvar = 0;
  do_scores = 0;
  do_CI_load = 0;
  do_CI_percvar = 0;

  if (nargout > 1)
    do_percvar = 1;
  end;
  if (nargout > 2)
    do_scores = 1;
  end;
  if (nargout > 3)
    do_CI_load = 1;
  end;
  if (nargout > 4)
    do_CI_percvar = 1;
  end;

  load = zeros(iter,nvars*ndf);   % Allocate loadings matrix
  if (do_percvar)
    percvar = zeros(iter,ndf);      % Allocate percent-variance matrix
  end;

  for it = 1:iter                 % Randomized discriminant analysis
    vincl = 1:nvars;              % List of vars to include
    for i = 1:(nvars-nv(n))       % Random deletions
      r = ceil(rand*(nvars-i+1));
      vincl(r) = [];
    end;
    Xi = X(:,vincl);              % Reduced data matrix

    mean_W = (G'*G)\G'*Xi;           % Within-group means
    devs_W = Xi - G*mean_W;          % Within-group deviations
    W = devs_W'*devs_W;              % Within-group SSCP matrix

    devs_T = Xi - ones(nobs,1)*mean(Xi);  % Total deviations
    T = devs_T'*devs_T;              % Total SSCP matrix

    B = T - W;                       % Among-sample SSCP matrix
    [evects,evals] = geneigen(B,W);  % Generalized eigen analysis
    evects = real(evects);           % If complex, convert to real
    evals = real(evals);

    if (do_percvar)
      percvari = 100*evals/sum(evals);  % Percents variance
      percvari = percvar(1:ndf,1);      % Retain subset
    end;

    scoresi = score(Xi,evects,ndf);   % Scores for subset of DFs
    scoresi = zscore(scoresi);         % Standardize scores by DF
    loadi = corr(Xi,scoresi);          % Recalc loadings as correlations with scores

    if (it>1)                        % After first iteration
      for d = 1:ndf                  %   check results against first
        if (corr(loadi(:,d),init_load(:,d))<0)
          loadi(:,d) = -loadi(:,d);
        end;
      end;
    else
      init_load = loadi;
    end;

    load(it,:) = loadi(:)';
    if (do_percvar)
      percvar(it,:) = percvari;
    end;
  end;



  return;
