% DISCRIMB: Objective function for bootstrapping discrim().  Returns a single 
%           row vector containing loadings and percvar values.
%
%     Usage: retstr = discrimb(X,grps,not_used1,not_used2,...
%                              ndf,loadtype,norescale,origloadings)
%
%         X =             [n x p] data matrix (obs x vars).
%         grps =          row or column vector of group identifiers.
%         ndf =           number of leading discriminant functions for
%                           which scores are desired (default = groups-1).
%         loadtype =      optional boolean flag indicating the scaling for the 
%                           loadings: 
%                             0: vector correlations [default];
%                             1: regression)coefficients;
%                             2: squared loadings sum to unity.
%         origloadings =  [p x ndf] matrix of loadings from original analysis.
%         --------------------------------------------------------------------
%         retstr =        row vector containing loadings and percvar results.
%

% RE Strauss, 6/3/99
%   4/29/00 - correct use of 'loadtype' flag.
%   5/4/00 -  call loadscrs() for loadings and scores.

function retstr = discrimb(X,grps,nu1,nu2,ndf,loadtype,norescale,origloadings)
  nobs = length(grps);

  G = design(grps);                 % ANOVA-type design matrix
  mean_W = (G'*G)\G'*X;             % Within-group means
  devs_W = X - G*mean_W;            % Within-group deviations
  W = devs_W'*devs_W;               % Within-group SSCP matrix

  devs_T = X - ones(nobs,1)*mean(X);  % Total deviations
  T = devs_T'*devs_T;               % Total SSCP matrix
  B = T - W;                        % Among-sample SSCP matrix

  [evects,evals] = geneigen(B,W);   % Generalized eigen analysis
  evects = real(evects);            % If complex, convert to real
  evals =  real(evals);
  evals = max([evals'; zeros(size(evals'))])';  % Ignore negative eigenvalues

  percvar = 100*evals/sum(evals);   % Percents variance
  percvar = percvar(1:ndf,1);       % Retain subset

  if (norescale)                    % Loadings and scores
    [loadings,scores] = loadscrs(X,evects,ndf,loadtype);
  else
    [loadings,scores] = loadscrs(X,evects,ndf,loadtype,grps);  
  end;

  for d = 1:ndf                     % Check if direction consistent with original
    if (corr(loadings(:,d),origloadings(:,d),2)<0)
      loadings(:,d) = -loadings(:,d);
    end;
  end;

  retstr = [loadings(:)' percvar'];

  return;

