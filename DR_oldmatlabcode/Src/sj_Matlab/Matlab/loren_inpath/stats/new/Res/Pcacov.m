% PCACOV: Returns the full set of sorted eigenvectors and eigenvalues, and
%          a subset of PCA scores, from the covariance matrix of a data
%          matrix X.  Optionally bootstraps the loadings and eigenvalues.
%
%     Syntax: [loadings,percvar,scores,fscores,CI_loadings,CI_percvar]
%              = pcacov(X,{npc},{loadtype},{Xf},{iter},{CI_level},{nowarn})
%
%         X =           [n x p] data matrix.
%         npc =         optional number of leading principal components to be 
%                         returned [default = number of significant eigenvalues 
%                         based on broken-stick null model].
%         loadtype =    optional boolean flag indicating the scaling for the 
%                         loadings: 
%                           0: vector correlations [default];
%                           1: regression coefficients;
%                           2: squared loadings sum to unity.
%         Xf =          optional [m x p] matrix of observations to be "floated" 
%                         onto the principal components.
%         iter =        optional number of bootstrap iterations for confidence 
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

% RE Strauss, 11/21/95
%   9/14/99 -   misc changes for Matlab v5.
%   11/21/99 -  updated documentation; 
%                 added capability of floating obs on PCs;
%                 implemented BOOTSTRP for bootstrapping.
%   1/15/00 -   allow for fewer observations than variables by estimating 
%                 loadings as regression coefficients.
%   5/2/00 -    isolated scores & loadings in LOADSCRS.
%   6/12/00 -   added error message for missing data.
%   12/11/01 -  added flag to supress warning messages.

function [loadings,percvar,scores,fscores,CI_loadings,CI_percvar] ...
            = pcacov(X,npc,loadtype,Xf,iter,CI_level,nowarn)

  if (nargin < 2) npc = []; end;
  if (nargin < 3) loadtype = []; end;
  if (nargin < 4) Xf = []; end;
  if (nargin < 5) iter = []; end;
  if (nargin < 6) CI_level = []; end;
  if (nargin < 7) nowarn = []; end;

  default_loadtype = 0;
  default_iter = 0;
  default_CI_level = 95;
  default_nowarn = 0;

  if (isempty(loadtype))
    loadtype = default_loadtype;
  end;
  if (isempty(iter))
    iter = default_iter;
  end;
  if (isempty(CI_level))
    CI_level = default_CI_level;
  end;
  if (isempty(nowarn))
    nowarn = default_nowarn;
  end;

  if (CI_level > 1)
    CI_level = CI_level/100;
  end;
  alpha = 1-CI_level;

  if (misscheck(X,Xf))
    error('  PCACOV: input matrix contains missing values');
  end;

  [N,P] = size(X);
  if (N < 3)
    error('  PCACOV: too few observations.');
  end;
  if (P < 1)
    error('  PCACOV: too few variables.');
  end;
  if (N < P & ~nowarn)
    if (iter)
      disp('  PCACOV warning: fewer observations than variables;');
      disp('                  bootstrapping not performed.');
    else
      disp('  PCACOV warning: fewer observations than variables.');
    end;
  end;

  if (~isempty(npc))
    if (npc > P)
      npc = P;
    end;
  else
    npc = P;
  end;

  if (~isempty(Xf))
    [M,Pf] = size(Xf);
    if (Pf ~= P)
      error('  PCACOV: floated data matrix must have same vars as data matrix.');
    end;
  end;

  % PCA
  covmat = cov(X);                        % Covariance matrix

  reduced_rank = 0;
  if (N >= P)                             % PCA
    [evects,evals] = eigen(covmat);
  else                                    % If too few obs,   
    v = sort(steprank(covmat));             % Find best subset of variables
    lenv = length(v);
    [ev,evals] = eigen(covmat(v,v));        % PCA of subset of vars
    s = score(X(:,v),ev);                   % Scores on subset of vars
    b = linregr(s,X);                       % Regress data on scores
    evects = b(2:lenv+1,:)';                % Evects are regression coeffs
    reduced_rank = 1;

    if (lenv < npc)
      if (~nowarn)
        disp('  PCACOV warning: due to singular covariance matrix, only');
        disp(['    ', sprintf('%d ',lenv), 'components can be returned.']);
      end;
      npc = lenv;
    end;
  end;

  if (isempty(npc))
    [pe,npc] = brokestk(P,sum(evals),evals);
    if (npc==0)
      if (~nowarn)
        disp('  PCACOV warning: no eigenvalues reported as significant.');
      end;
      npc = 1;
    end;
  end;

  percvar = 100 * evals / sum(evals);     % Percents variance
  percvar = percvar(1:npc);               % Retain subset

  [loadings,scores] = loadscrs(X,evects,npc,loadtype);

  if (~isempty(Xf))                       % Float additional obs onto PCs
    fscores = score(Xf,evects,npc);
  else
    fscores = [];
  end;

  % Bootstrap PCA

  CI_loadings = [];
  CI_percvar = [];

  if (iter & ~reduced_rank)
    ci = bootstrp('pcacovb',1,iter,alpha,X,[],0,npc,loadtype,loadings);

    nn = P*npc;                           % Reshape into CI matrices
    CI_loadings = [reshape(ci(1,1:nn)',P,npc) ...
                   reshape(ci(2,1:nn)',P,npc)];

    c = [];
    for i = 1:npc
      c = [c i:npc:2*npc];
    end;
    CI_loadings = CI_loadings(:,c);

    CI_percvar =  [reshape(ci(1,nn+1:nn+npc),npc,1) ...
                   reshape(ci(2,nn+1:nn+npc),npc,1)];
  end;

  return;

