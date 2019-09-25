% DISCRIMF: Calculating function for discrim().
%
%     Usage: [loadings,percvar,scores,fscores,B,W] = ...
%                                     discrimf(X,grps,Xf,ndf,loadtype,norescale)
%
%         X =         [n x p] data matrix (obs x vars).
%         grps =      row or column vector of group identifiers.
%         Xf =        optional [m x p] matrix of observations to be "floated" 
%                       onto the discriminant functions.
%         ndf =       number of leading discriminant functions for
%                       which scores are desired (default = groups-1).
%         loadtype =  boolean flag indicating the scaling for the 
%                       loadings: 
%                         0: vector correlations [default];
%                         1: regression)coefficients;
%                         2: squared loadings sum to unity.
%         norescale = boolean flag indicating that scores are not to 
%                       be rescaled [default = 0].
%         ----------------------------------------------------------------------
%         loadings =  [p x ndf] matrix of discriminant-function
%                       loadings (columns) as vector correlations or 
%                       regression coefficients.
%         percvar =   column vector of percents of total variance explained
%                       for discriminant functions.
%         scores =    [p x ndf] matrix of discriminant scores (columns).
%         fscores -   [m x npc] matrix of discriminant scores for "floated" obs.
%         B =         among-sample SSCP matrix.
%         W =         within-group SSCP matrix.
%

% RE Strauss, 6/3/99
%   11/29/99 - changed calling sequence.
%   4/29/00 -  correct usage of 'loadtype' flag.
%   5/4/00 -   call loadscrs() for loadings and scores; added 'norescale' option.
%   11/14/00 - correct problem with rescaling of floated scores.

function [loadings,percvar,scores,fscores,B,W] = ...
                      discrimf(X,grps,Xf,ndf,loadtype,norescale)
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

  if (norescale)
    [loadings,scores] = loadscrs(X,evects,ndf,loadtype);
  else
    [loadings,scores,meanscr,wvar] = loadscrs(X,evects,ndf,loadtype,grps);
  end;

  if (~isempty(Xf))                 % Scores to be 'floated'
    fscores = score(Xf,evects,ndf);
    nf = size(Xf,1);
    if (~norescale)
      fscores = (fscores-ones(nf,1)*meanscr)./(ones(nf,1)*sqrt(wvar));
    end;
  else
    fscores = [];
  end;

  for d = 1:ndf                     % Check if DF direction positive
    if (sum(loadings(:,d))<0)
      loadings(:,d) = -loadings(:,d);
      scores(:,d) = -scores(:,d);
      if (~isempty(fscores))
        fscores(:,d) = -fscores(:,d);
      end;
    end;
  end;

  return;

