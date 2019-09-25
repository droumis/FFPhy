% PCACORR: Returns the full set of sorted eigenvectors and eigenvalues, and
%          a subset of PCA scores, from the correlation matrix of a data
%          matrix X.  Optionally bootstraps the loadings and eigenvalues.
%
%     Syntax: [loadings,percvar,scores,CI_loadings,CI_percvar]
%              = pcacorr(X,{npc},{loadtype},{Xf},{boot_iter},{CI_level})
%
%         X =           [n x p] data matrix.
%         npc =         optional number of leading principal components to be 
%                         returned [default=3].
%         loadtype =    optional boolean flag indicating the scaling for the 
%                         loadings: 
%                           0: vector correlations [default];
%                           1: eigenvector (regression) coefficients;
%                           2: squared loadings sum to unity.
%         Xf =          optional [m x p] matrix of observations to be "floated" 
%                         onto the principal components.
%         boot_iter =   optional number of bootstrap iterations for confidence 
%                         intervals [default=0].
%         CI_level =    optional width of confidence intervals [default=95].
%         nowarn =      optional boolean flag indicating, if true, that warning
%                         messages are to be suppressed [default=0].
%         ---------------------------------------------------------------------
%         loadings =    [p x npc] matrix of principal components (columns).
%         percvar =     [p x 1] vector of percentages of variance accounted for
%                         by principal components.
%         scores =      [n x npc] matrix of unscaled PCA scores (columns).
%         fscores -     [m x npc] matrix of PCA scores for "floated" obs.
%         CI_loadings = [p x 2*npc] matrix of CI% asymmetric confidence limits
%                         of loadings, two columns per component (low,high).
%         CI_percvar =  [p x 2] matrix of CI% asymmetric confidence limits
%                         of percents variance-explained.
%

% RE Strauss, 11/10/97
%   9/19/99 - updated handling of default input arguments.
%   11/21/99 - altered to call PCACOV with standardized variables.
%   12/11/01 -  added flag to supress warning messages.

function [loadings,percvar,scores,fscores,CI_loadings,CI_percvar] ...
            = pcacorr(X,npc,loadtype,Xf,iter,CI_level,nowarn)

  if (nargin < 2) npc = []; end;
  if (nargin < 3) loadtype = []; end;
  if (nargin < 4) Xf = []; end;
  if (nargin < 5) iter = []; end;
  if (nargin < 6) CI_level = []; end;
  if (nargin < 7) nowarn = []; end;

  X = zscore(X);                        % Standardize variables
  if (~isempty(Xf))
    Xf = zscore(Xf);
  end;

  if (misscheck(X,Xf))
    error('  PCACORR: input matrix contains missing values');
  end;

  [loadings,percvar,scores,fscores,CI_loadings,CI_percvar] ...
            = pcacov(X,npc,loadtype,Xf,iter,CI_level,nowarn);

  return;

