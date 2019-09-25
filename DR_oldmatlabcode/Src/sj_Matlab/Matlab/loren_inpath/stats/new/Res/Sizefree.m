% SIZEFREE: Size-invariant discriminant analysis
%
%     Syntax: [loadings,percvar,scores,fscores,D2,R,wload,wperc,wsize,S]
%              = sizefree(X,grps,{ndf},{loadtype},{kindsize},{kindregr},{Xf},{iter},{CI_level})
%
%        X =           [n x p] data matrix (obs x vars).
%        grps =        row or column vector of group identifiers.
%        ndf =         optional number of leading discriminant functions for
%                        which scores are desired [default=groups-1, max 10].
%        loadtype =    optional boolean flag indicating the scaling for the 
%                        loadings: 
%                          0: vector correlations [default];
%                          1: regression)coefficients;
%                          2: squared loadings sum to unity.
%        kindsize =    kind of size vector: 'w' for within-group or
%                        or 'a' for among-group [default = 'w'].
%        kindregr =    regression model used: 0 = major axis [default], 1 = predictive.
%        Xf =          optional [m x p] matrix of observations to be "floated" 
%                        onto the discriminant functions.
%        iter =        number of bootstrap iterations [default=0].
%        CI_level =    percent width of confidence intervals [default=95].
%        ------------------------------------------------------------------------------
%        loadings =    [p x ndf] matrix of discriminant-function
%                        coefficients (columns) as vector correlations
%                        ([p x 3*ndf] if bootstrapped).
%        percvar =     [ndf x 1] vector of percents of total variance explained
%                        for discriminant functions ([ndf x 3] if bootstrapped).
%        scores =      [n x ndf] matrix of discriminant-function scores
%                        (columns)  (can't be bootstrapped).
%        fscores =     [m x ndf] matrix of "floated" observations.
%        D2 =          [g x g] symmetric Mahalanobis distance matrix of residuals 
%                        ([g x 3*g] if bootstrapped).
%        R =           [n x p] matrix of size-invariant residuals (can't 
%                        be bootstrapped).
%        wload =       [p x 1] vector of size-vector loadings ([p x 3] if 
%                        bootstrapped).
%        wperc =       percent total variance for size vector ([1 x 3] if 
%                        bootstrapped).
%        wsize =       [n x 1] vector of within-group size scores (can't 
%                        be bootstrapped).
%        S =           [n x 1] vector of among-group size scores (can't 
%                        be bootstrapped).
%
%        Note: bootstrapping currently provides only confidence limits, 
%              not significance levels or power estimates, in the form:
%
%              output_matrix = [point estimates, lower CL, upper CL]
%

% RE Strauss, 6/12/95
%   7/17/98 - added optional major-axis regression for estimating residuals.
%  11/29/99 - reversed X and grps in calling sequence.
%   6/13/00 - added check for missing data.

function [loadings,percvar,scores,fscores,D2,R,wload,wperc,wsize,S] ...
              = sizefree(X,grps,ndf,loadtype,kindsize,kindregr,Xf,iter,CI_level)

  if (nargin < 3) ndf = []; end;
  if (nargin < 4) loadtype = []; end;
  if (nargin < 5) kindsize = []; end;
  if (nargin < 6) kindregr = []; end;
  if (nargin < 7) Xf = []; end;
  if (nargin < 8) iter = []; end;
  if (nargin < 9) CI_level = []; end;

  [nobs,nvars] = size(X);               % Numbers of observations & variables
  [nfobs,nfvars] = size(Xf);            % Numbers of floated obs & vars
  index = uniquef(grps);
  ngrps = length(index);                % Number of groups

  if (~isempty(Xf) & nvars~=nfvars)
    error('  SIZEFREE: numbers of variables for obs and floated obs must be identical.');
  end;

  if (misscheck(X,grps,Xf))
    error('  SIZEFREE: data matrix or grouping vector contains missing data.');
  end;

  default_ndf = ngrps-1;                % Defaults for input arguments
  default_iter = 0;
  default_CI_level  = 95;
  default_kindsize = 'w';               % = within
  default_loadtype = 0;                 % = vector correlations
  default_kindregr = 0;                 % = major-axis regression

  if (isempty(ndf))                     % Default input arguments
     ndf = min(default_ndf,10);
  end;
  if (isempty(kindsize))
    kindsize = default_kindsize;
  end;
  if (kindsize == 'w')
    within = 1;
  else
    within = 0;
  end;
  if (isempty(kindregr))
    kindregr = default_kindregr;
  end;
  if (isempty(iter))
    iter = default_iter;
  end;
  if (isempty(CI_level))
    CI_level = default_CI_level;
  end;


  % Calculate single solution

%        outmat =      boolean vector indicating results to be returned:
%                         1) loadings: (vector correlations)        [p x ndf] 
%                         2) percvar:  percents total variance      [ndf x 1] 
%                         3) scores:   DF scores                    [n x ndf]
%                         4) fscores:  floated DF scores            [m x ndf]
%                         5) D2:       Mahal distances of residuals [g x g]
%                         6) R:        size-invariant residuals     [n x p]
%                         7) wload:    size-vector loadings         [p x 1]
%                         8) wperc:    percent size-vector variance [1 x 1]
%                         9) wsize:    within-group size scores     [n x 1]
%                        10) S:        among-group size scores      [n x 1]

  outmat = [1 zeros(1,9)];
  outmat(2:nargout) = ones(1,nargout-1);

  outsize = [nvars*ndf, ndf, nobs*ndf, nfobs*ndf, ngrps^2, nobs*nvars, nvars, nobs, nobs];
  i = find(~outmat);
  if (~isempty(i))
    outsize(i) = zeros(1,length(i));
  end;
  i = find(~outsize);
  if (~isempty(i))
    outmat(i) = zeros(1,length(i));
  end;

  solution = sizefref(X,grps,Xf,[],[],ngrps,ndf,within,outmat,outsize,loadtype,kindregr);

  [loadings,percvar,scores,fscores,D2,R,wload,wperc,wsize,S] ...
                                = sizefrep(solution,outmat,outsize,nobs,nvars,ngrps,ndf);


  % Bootstrapping

  if (iter > 0)
    outmat([3 5 8 9]) = zeros(1,4);   % Omit scores from bootstrapping
    if (ndf==1)                       % Omit percvar if = 100%
      outmat(2) = 0;
    end;
    outsize(~outmat) = zeros(1,sum(~outmat));

    ci = bootstrp('sizefref',[1 0 0],iter,1-CI_level,X,grps,Xf,0, ...
                  ngrps,ndf,within,outmat,outsize,loadtype,kindregr);
    [loadings1,percvar1,scores1,D21,R1,wload1,wperc1,wsize1,S1] ...
               = sizefrep(ci(1,:),outmat,nobs,nvars,ngrps,ndf);
    [loadings2,percvar2,scores2,D22,R2,wload2,wperc2,wsize2,S2] ...
               = sizefrep(ci(2,:),outmat,nobs,nvars,ngrps,ndf);

    if (outmat(1))
      loadings = [loadings loadings1 loadings2];
    end;
    if (outmat(2))
      percvar = [percvar percvar1 percvar2];
    end;
    if (outmat(4))
      D2 = [D2 D21 D22];
    end;
    if (outmat(6))
      wload = [wload wload1 wload2];
    end;
    if (outmat(7))
      wperc = [wperc wperc1 wperc2];
    end;
  end;

  return;

