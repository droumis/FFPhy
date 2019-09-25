% VECTCORR: Calculates the vector correlations among a set of vectors (given 
%       by coefficients), or between a reference vector and a set of vectors,
%       or among all possible pairs of two sets of vectors.
%       Returns a single correlation for two vectors.
%          Missing coefficents are replaced by means, by vector.
%          Provides 2-tailed probabilities, either as asymptotic small-sample
%       estimates or via independent randomization of variables.  Probabilities 
%       are NOT simultaneous (see CORRPROB).
%          Provides confidence intervals for correlation coefficients, either as
%       asymptotic small-sample estimates or via bootstrapping original
%       observations.  Confidence intervals are NOT simultaneous (see CORRCI).
%
%     Syntax: [r,prob,CI_low,CI_high] = vectcorr(X1,{X2},{iter},{CI_level})
%
%       X1 =        [N x P1] data matrix.  If P1=1, it is taken to be the
%                     reference vector.
%       X2 =        optional [N x P2] data matrix.  If P2=1, it is taken to be
%                     the reference vector.  If X2 is a scalar, it is assumed 
%                     that X2 has been omitted and special, iter and CI_level 
%                     are in argument positions 2-4.
%       iter =      optional number of iterations for randomized probabilities  
%                     and confidence intervals [default = 0].
%       CI_level =  optional percentage width of confidence intervals 
%                     [default = 95].
%       -------------------------------------------------------------------------
%       r =         [P1 x P1] symmetric correlation (if only X1 is specified);
%                   [P1 x 1] vector of correlations (if X1 & ref are specified);
%                   [1 x P2] vector of correlations (if ref & X2 are specified);
%                   [P1 x P2] rectangular matrix of correlations (if X1 & X2 are
%                     both specified);
%                   [1 x 1] scalar correlation (if X1 is an [N x 2] matrix or
%                     ref1 & ref2 are both specified).
%       prob =      corresponding probabilities, either asymptotic (small-sample)
%                     if iter=0, or via randomization if iter>0.
%       CI_low =    corresponding lower confidence limits, either asymptotic
%                     (small-sample) if iter=0, or via bootstrapping
%                     if iter>0.
%       CI_high =   corresponding upper confidence limits.
%

% Note: Any changes to the documentation of this function should be made in parallel
%       to CORR() and RANKCORR().

% RE Strauss
%   9/9/96 - missing values replaced by means for each variable.
%   5/25/99 - modified to call corr()

function [r,prob,CI_low,CI_high] = vectcorr(X1,X2,iter,CI_level)
  if (nargin < 2) X2 = []; end;
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  special = 2;
  if (nargout < 2)
    r = corr(X1,X2,special);
  elseif (nargout < 3)
    [r,prob] = corr(X1,X2,special,iter);
  else
    [r,prob,CI_low,CI_high] = corr(X1,X2,special,iter,CI_level);
  end;

  return;
