% MORPHINTI:  Calculates the value of the morphological-integration index from 
%             the correlation or covariance matrix of a data matrix.  
%             Uses Leamy's (1994) modification of a measure of integration 
%             proposed by Cheverud, Rutledge & Atchley (1983:897) (the adjusted 
%             variance among the eigenvalues).  Eigenvalues are scaled to a sum 
%             of one before integration is estimated, consistent with sets of 
%             eigenvalues estimated from correlation matrices.  'Size-free' 
%             patterns of integration are assessed by ignoring the first 
%             eigenvalue, under the assumption that all size variation is 
%             described by PC1, and scaling the remaining eigenvalues to a mean 
%             of one.  The index varies from zero (no correlations among 
%             characters) to unity (singular covariance or correlation matrix).
%               Called by MORPHINT.
%
%     Usage:  [I,evals,pc1] = morphinti(X,sizefree,usecorr)
%
%           X =         [n x p] data matrix for single group of specimens.
%           sizefree =  boolean variable indicating whether (=1) or not (=0) to 
%                         ignore PC1 in estimation of integration indices 
%                         [default = 0].
%           usecorr =   boolean variable indicating (=1) that principal 
%                         component analyses are to done on the correlation 
%                         matrices of characters [default = covariance matrices].
%           ---------------------------------------------------------------------
%           I =         integration index.
%           evals =     column vector of eigenvalues (including first, even if I 
%                         is sizefree).
%           pc1 =       PC1 loadings (even if I is sizefree).
%

% RE Strauss, 6/2/01 (extracted from morphint.m)
%   7/12/01 - export eigenvalues and PC1 loadings.

function [I,evals,pc1] = morphinti(X,sizefree,usecorr)
  if (nargin < 2) sizefree = []; end;
  if (nargin < 3) usecorr = []; end;

  getpc1 = 0;
  if (nargout > 2)
    getpc1 = 1;
  end;

  if (isempty(sizefree))                  % Default argument values
    sizefree = 0;
  end;
  if (isempty(usecorr))
    usecorr = 0;
  end;

  [nobs,nvars] = size(X);

  if (usecorr)                            % Correlation or covariance matrix
    c = corrcoef(X);
  else
    c = cov(X);
  end;

  [evects,evals] = eigen(c);              % Eigenanalysis
  e = evals;

  pc1 = [];
  if (getpc1)
    pc1 = loadscrs(X,evects,1,0);
  end;

  if (sizefree)                           % If sizefree, ignore max eigenvalue
    e = e(2:end);
  end;
  e = e./mean(e);                         % Standardize to mean=1

  I = sqrt(var(e)/nvars);                 % Integration coefficient
% I = sqrt(var(e)/sum(e));

  return;

