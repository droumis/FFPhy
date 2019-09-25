% LOADSCRS: Calculate loadings and scores.
%
%     Usage: [loadings,scores,meanscr,wvar] = loadscrs(X,evects,nvect,loadtype,{zs})
%
%         X =         [n x p] data matrix.
%         evects =    [p x nvect] matrix of eigenvectors.
%         nvect =       optional number of leading eigenvectors to be 
%                       returned.
%         loadtype =  optional boolean flag indicating the scaling for the 
%                       loadings: 
%                         0: vector correlations [default];
%                         1: regression coefficients;
%                         2: squared loadings sum to unity.
%         zs =        optional boolean flag indicating, if true, that scores 
%                       are to be standardized to zero mean and unit variance 
%                       [default = 0].  
%                     If loadtype = 1, this is done before regression 
%                       coefficients are estimated.
%                     If zs is a group-membership vector rather than a scalar, 
%                       scores are standardized to zero mean and within-group 
%                       variance.
%         ---------------------------------------------------------------------
%         loadings =  [p x nvect] matrix of loadings.
%         scores =    [n x nvect] matrix of scores.
%         meanscr =   vector of mean scores before scaling (if zs>0).
%         wvar =      vector of within-group variances in scores before scaling 
%                       (if zs>0).
%

% RE Strauss, 5/2/00
%   5/3/00 - added error for invalid 'loadtype' flag.
%   5/4/00 - added option to standardize scores; return meanscr and wvar.

function [loadings,scores,meanscr,wvar] = loadscrs(X,evects,nvect,loadtype,zs)
  if (nargin < 5) zs = []; end;

  if (isempty(zs))
    zs = 0;
  end;
  if (isvector(zs))
    grps = zs;
    zs = 2;
  end;

  [n,p] = size(X);
  scores = score(X,evects,nvect);           % Scores for subset of PCs

  meanscr = [];
  wvar = [];

  switch(zs)                                % Standardize scores
    case 0
    case 1,
      scores = zscore(scores);

    case 2,
      Cp = covpool(scores,grps,1);
      meanscr = mean(scores);                   % Mean scores before scaling
      wvar = diag(Cp)';                         % Pooled within-group variances
      scores = (scores-ones(n,1)*meanscr)./(ones(n,1)*sqrt(wvar));

    otherwise
      error('  LOADSCRS: invalid zs flag');
  end;

  switch(loadtype)                          % Loadings
    case 0,                                   % Vector correlations
      loadings = corr(X,scores);            

    case 1,                                   % Regression coefficients
      loadings = zeros(p,nvect);
      for ipc = 1:nvect
        S = [ones(n,1) scores(:,ipc)];
        for ip = 1:p
          y = X(:,ip);
          b = inv(S'*S)*S'*y;
          loadings(ip,ipc) = b(2);
        end;
      end;

    case 2,
      loadings = sumsqscale(evects(:,1:nvect)); % Squared coeffs sum to unity

    otherwise
      error('  LOADSCRS: invalid loadings type');
  end;

  return;
