% PCACONSTR: Constrained principal component analysis of a residual covariance 
%            matrix of Y, after sweeping out a set of independent variables X.
%            Removes observations having missing values.
%              Example of use: PCA of morphological variables, independent of a
%            set of ecological variables.
%
%     Usage: [loadings,percvar,scores,b,stats,CI_loadings,CI_percvar] = ...
%                                       pcaconstr(X,Y,{npc},{loadtype},{iter})
%
%         X =           [n x q] matrix of independent variables.
%         Y =           [n x p] matrix of dependent variables.
%         npc =         optional number of leading principal components to be 
%                         returned [default = all].
%         loadtype =    optional boolean flag indicating the scaling for the 
%                         loadings: 
%                           0: vector correlations [default];
%                           1: regression coefficients;
%                           2: squared loadings sum to unity.
%         iter =        optional number of bootstrap iterations [default=0].
%         --------------------------------------------------------------------
%         loadings =    [p x npc] matrix of principal components (columns).
%         percvar =     [p x 1] vector of percents variance-explained
%                         for principal components.
%         scores =      [n x npc] matrix of PCA scores (columns).
%         b =           [2 x p] matrix of parameter estimates from multiple 
%                         regression.
%         stats =       [6 x p] matrix of regression statistics from multiple 
%                         regression [see linregr()].
%         CI_loadings = [p x 2*npc] matrix of CI% asymmetric confidence limits
%                         of loadings, two columns per component (low,high).
%         CI_percvar =  [p x 2] matrix of CI% asymmetric confidence limits
%                         of percents variance-explained.
%

% Little, RJA & DB Rubin. 1987. Statistical Analysis with Missing Data. Wiley.
%   Section 6.5, pp. 112-115.

% RE Strauss, 6/1/00
%   11/28/00 - changed default npc to all rather than 3.

function [loadings,percvar,scores,b,stats,CI_loadings,CI_percvar] = ...
            pcaconstr(X,Y,npc,loadtype,iter)

  if (nargin < 3) npc = []; end;
  if (nargin < 4) loadtype = []; end;
  if (nargin < 5) iter = []; end;

  rx = rowsum(X);                         % Remove obs having missing data
  ry = rowsum(Y);
  if (any(~isfinite(rx)) | any(~isfinite(ry)))
    i = find(isfinite(rowsum([rx,ry])));
    X = X(i,:);
    Y = Y(i,:);
  end;

  [n,q] = size(X);
  [ny,p] = size(Y);

  if (n~=ny)
    error('  PCACONSTR: input matrices must have same number of rows.');
  end;

  if (isempty(npc))
    npc = p;
  end;
  if (isempty(loadtype))
    loadtype = 0;
  end;
  if (isempty(iter))
    iter = 0;
  end;

  iv = 1:size(X,2);                       % Indices of independent variables

  [soln,scores,b,stats] = pcaconstrf([X Y],[],[],[],iv,npc,loadtype);
  loadings = reshape(soln(1:(p*npc))',p,npc);
  percvar = soln((p*npc+1):length(soln))';

  if (iter)

  else
    CI_loadings = [];
    CI_percvar = [];
  end;

  return;

