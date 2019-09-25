% DISCRIM: Canonical variate (discriminant) analysis.
%
%     Usage: [loadings,percvar,scores,fscores,CI_loadings,CI_percvar]
%              = discrim(X,grps,{ndf},{loadtype},{Xf},{norescale},{iter},{CI_level})
%
%         X =           [n x p] data matrix (obs x vars).
%         grps =        row or column vector of group identifiers.
%         ndf =         optional number of leading discriminant functions for
%                         which scores are desired (default = groups-1).
%         loadtype =    optional boolean flag indicating the scaling for the 
%                         loadings: 
%                           0: vector correlations [default];
%                           1: regression coefficients;
%                           2: squared loadings sum to unity.
%         Xf =          optional [m x p] matrix of observations to be "floated" 
%                         onto the discriminant functions.
%         norescale =   optional boolean flag indicating that scores are not to 
%                         be rescaled [default = 0].
%         iter =        optional number of bootstrap iterations [default=0].
%         CI_level =    optional percent width of confidence intervals [default=95].
%         --------------------------------------------------------------------------
%         loadings =    [p x ndf] matrix of discriminant-function
%                         loadings (columns) as vector correlations or 
%                         regression coefficients.
%         percvar =     column vector of percents of total variance explained
%                         for discriminant functions.
%         scores =      [p x ndf] matrix of discriminant scores (columns).
%         fscores -     [m x npc] matrix of discriminant scores for "floated" obs.
%         CI_loadings = [p x 2*ndf] matrix of CI% confidence limits
%                         (asymmetric) of vector-correlation loadings, two columns
%                         per discriminant function (low, high).
%         CI_percvar =  [p x 2] matrix of CI% confidence limits (asymmetric)
%                         of percents variance-explained.
%

% RE Strauss, 6/5/95
%   11/9/98 -   warning message for singular W or B matrices.
%   6/3/99 -    major rewrite.
%   11/21/99 -  corrected sequence of CI_loadings cols.
%   11/25/99 -  added option for 'floated' scores.
%   11/29/99 -  reversed X and grps in calling sequence; also in discrimf().
%   6/13/00 -   added check for missing data.
%   4/6/02 -    remove error checks for reduced-rank matrices;
%               set number of functions returned as the minimum of the requested
%               ndf, k-1, and p.
%   6/25/02 -   undid previous change: returns 'ndf' functions.

function [loadings,percvar,scores,fscores,CI_loadings,CI_percvar] ...
            = discrim(X,grps,ndf,loadtype,Xf,norescale,iter,CI_level)

  if (nargin < 3) ndf = []; end;
  if (nargin < 4) loadtype = []; end;
  if (nargin < 5) Xf = []; end;
  if (nargin < 6) norescale = []; end;
  if (nargin < 7) iter = []; end;
  if (nargin < 8) CI_level = []; end;

  CI_percvar = [];                % Allocate optional return arguments
  CI_loadings = [];

  ngrps = length(unique(grps));   % Number of groups
  [nobs,nvars] = size(X);         % Numbers of observations & variables

  if (length(grps) ~= nobs)
    error('  DISCRIM: Group vector and data matrix not compatible');
  end;
  
  if (misscheck(X,grps,Xf))
    error('  DISCRIM: data matrix or grouping vector contains missing data.');
  end;

  if (isempty(ndf))
    ndf = min([ngrps-1,nvars]);
  end;
  if (isempty(loadtype))
    loadtype = 0;
  end;
  if (isempty(norescale))
    norescale = 0;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  end;

  if (CI_level > 1)
    CI_level = CI_level/100;
  end;
  alpha = 1-CI_level;
  
%   ndf = min([ndf,ngrps-1,nvars]);         % Set number of functions returned

  % Discriminant analysis

  [loadings,percvar,scores,fscores,B,W] = discrimf(X,grps,Xf,ndf,loadtype,norescale);

  sizeW = size(W,1);
  rankW = rank(W);
  rankB = rank(B);

%   if (rankW < sizeW)
%     disp(sprintf('  DISCRIM warning: W matrix is singular (rank %1.0f < %1.0f).',...
%                     rankW, sizeW));
%   end;
%   if (rankB < ngrps-1)
%     disp(sprintf('  DISCRIM warning: B matrix is singular (rank %1.0f < %1.0f).',...
%                     rankW, ngrps-1));
%   end;

  % Bootstrap loadings and percvar

  if (iter)
    ci = bootstrp('discrimb',1,iter,alpha,X,grps,0,ndf,loadtype,norescale,loadings);

    nn = nvars*ndf;                   % Reshape into CI matrices
    CI_loadings = [reshape(ci(1,1:nn)',nvars,ndf) ...
               reshape(ci(2,1:nn)',nvars,ndf)];

    c = [];
    for i = 1:ndf
      c = [c i:ndf:2*ndf];
    end;
    CI_loadings = CI_loadings(:,c);

    CI_percvar = [reshape(ci(1,nn+1:nn+ndf),ndf,1) ...
                  reshape(ci(2,nn+1:nn+ndf),ndf,1)];
  end;
      
  return;
