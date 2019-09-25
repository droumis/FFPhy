% CANCORR:  Canonical correlation analysis.
%             Significance of canonical correlations is assessed by permuting 
%           the sequence of observations in one suite of variables, so that 
%           within-suite correlations are maintained but between-suite 
%           correlations are broken.  Because canonical correlations are always 
%           positive, the test is right-tailed.
%             Scores on canonical functions are standardized if loadtype==0.
%
%     Usage: [r,x_load,y_load,x_scores,y_scores,prr] 
%                              = cancorr(X,Y,{ncorr},{loadtype},{iter})
%
%         X =        [n x p] matrix for first suite of variables.
%         Y =        [n x q] matrix for second suite of variables.
%         ncorr =    number of canonical variates [default = min(p,q)].
%         loadtype = optional boolean flag indicating the scaling for the 
%                      loadings: 
%                        0: vector correlations [default];
%                        1: regression coefficients;
%                        2: scoring coefficients.
%         iter =     optional number of iterations for significance levels of
%                      canonical correlations [default = 0].
%         ------------------------------------------------------------------------
%         r =        [c x 1] vector of canonical correlations 
%                      [where c = min(p,q)].
%         x_load =   [p x ncorr] matrix of X variables (rows) with structural 
%                      correlations of canonical variates (cols) 
%                      for first set of variables.
%         y_load =   [q x ncorr] matrix of Y variables (rows) with structural 
%                      correlations of canonical variates (cols)
%                      for second set of variables.
%         x_scores = [n x ncorr] matrix of canonical-variate scores for first set 
%                      of variables.
%         y_scores = [n x ncorr] matrix of canonical-variate scores for second set 
%                      of variables.
%         prr =      [c x 1] vector of randomized significance levels of canonical 
%                      correlations.
%

% RE Strauss, 12/8/97
%   9/19/99 - updated handling of default input arguments.
%   11/6/00 - corrected problem with loadtype; conditions had been reversed;
%             adjust directions of vectors;
%             use switchem() function for switching matrix contents.
%   2/23/01 - correct problem with adjusting direction of vector for single 
%               variable.
%   11/7/01 - isolate computations in cancorrf();
%             add randomized significance testing of canonical correlations.

function [r,x_load,y_load,x_scr,y_scr,prr] = cancorr(X,Y,ncorr,loadtype,iter)
  if (nargin < 3) ncorr = []; end;
  if (nargin < 4) loadtype = []; end;
  if (nargin < 5) iter = []; end;

  [nobs,p] = size(X);
  [m,q] = size(Y);
  ncvar = min([q p]);                   % Number of pairs of canonical variables

  if (isempty(ncorr))
    ncorr = ncvar;
  end;
  if (isempty(loadtype))
    loadtype = 0;
  end;
  if (isempty(iter))
    iter = 0;
  end;

  if (q > p)                            % Enforce p >= q
    switch_sets = 1;  
  else
    switch_sets = 0;
  end;

  if (nobs ~= m)
    error('  CANCORR: input matrices must have same number of observations (=rows)');
  end;

  X = zscore(X);                        % Standardize variables
  Y = zscore(Y);

  if (switch_sets)                      % Switch sets of vars if p<q
    [X,Y] = switchem(X,Y);
    [p,q] = switchem(p,q);
  end;

  [r,x_load,y_load,x_scr,y_scr] = cancorrf(X,Y,ncorr,loadtype);

  prr = [];
  if (iter)                             % Significance levels of canonical corrs
    distrib = zeros(iter,ncorr);
    for it = 1:iter
      X = X(randperm(nobs),:);
      distrib(it,:) = cancorrf(X,Y,ncorr,loadtype)';
    end;
    prr = (sum(distrib >= ones(iter,1)*abs(r'))/iter)'; % Right-tailed
  end;

  if (switch_sets)                      % Re-switch sets of variables
    [x_scr,y_scr] = switchem(x_scr,y_scr);
    [x_load,y_load] = switchem(x_load,y_load);
    [p,q] = switchem(p,q);
  end;

  return;
