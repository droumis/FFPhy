% MORPHINTD:  Tests for a relationship between level of correlation and distance 
%             between two morphometric characters (e.g., measured between 
%             character midpoints) using Mantel's test.  Distances are assumed 
%             to represent the lower triangular portion of a distance matrix.  
%             Corrects pairwise distances by multiple regression on character 
%             magnitudes.  Because Mantel's procedure tests the equality of two 
%             distance matrices, correlations (r) are converted to 1-abs(r) as a 
%             measure of non-correlation.
%
%     Usage: [r,Z,pr,ncorr,pdist] = morphintd(X,dist,{iter},{noplot})
%
%             X =       [n x p] data matrix for n observations and p characters.
%             dist =    [n x p(1-p)/2] matrix of distances among characters for 
%                         each observation.
%             iter =    optional number of permutation iterations for Mantel's 
%                         test [default=0].
%             noplot =  optional boolean variable indicating, if true, that
%                         scatterplot is to be suppressed [default = 0].
%             ------------------------------------------------------------------
%             r =       correlation between matrix cells.
%             Z =       sum of cross products of matrix cells.
%             pr =      right-tailed significance level (p-value), either 
%                         asymptotic (if iter=0) or randomized.
%             ncorr =   non-correlation measures among characters.
%             pdist =   predicted distances among characters.
%

% RE Strauss, 6/2/01

function [r,Z,pr,corrinv,preddist] = morphintd(X,dist,iter,noplot)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) noplot = []; end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;

  [nobs,nchars] = size(X);
  [nobsd,ndists] = size(dist);

  if (nobs ~= nobsd)
    error('  MORPHINTD: number of observations must be same for both matrices.');
  end;
  if (ndists ~= nchars*(nchars-1)/2)
    error('  MORPHINTD: number of characters and distances not compatible.');
  end;

  corrs = corr(X);                      % Correlations among characters
  corrinv = 1-abs(corrs);               % Measure of non-correlation
  corrinv = putdiag(corrinv,0);         % Make diagonal elements exactly zero
  charmeans = means(X);                 % Character means across all individuals

  preddist = zeros(nchars,nchars);
  k = 0;
  for i = 2:nchars                      % Regress distances on character sizes
    for j = 1:i-1
      k = k+1;
      y = dist(:,k);
      x = X(:,[j,i]);
      im = find(isfinite(y) & isfinite(x(:,1)) & isfinite(x(:,2)));
      b = regress(y(im),x(im,:));
      preddist(i,j) = charmeans([j,i])*b;
    end;
  end;
  preddist = preddist + preddist';      % Reflect values
preddist = exp(preddist);
preddist = putdiag(preddist,0);

  [r,Z,pr] = mantel(corrinv,preddist,iter);   % Mantel's test

  if (~noplot)
    figure;
    tridist = trilow(preddist);
    tricorr = trilow(corrinv);
    plot(tridist,tricorr,'ko');
    putbnd(tridist,tricorr);
    putybnd(0,1);
    putxlab('Distance between characters');
    putylab('Non-correlation of characters');
  end;

  return;
