% KMEANGRP: Uses k-means clustering to determine the number of clusters.
%           Varies the numbers of clusters from 2 to a specified maximum and
%           evaluates the clustering solution based on the Krzanowski & Lai
%           (1988) sum-of-squares criterion.
%
%     Syntax: [optk,centr,clst,Cg] = kmeangrp(X,{maxk},{restarts},{doplot})
%
%           X =        [n x p] data matrix of coordinates. 
%           maxk =     optional maximum number of clusters fitted [default = 6].
%           restarts = optional number of restarts, after finding a new minimum sse,
%                        needed to end search for each value of k [default = 0].
%           doplot =   optional boolean flag producing plots if true [default = 0].
%           ------------------------------------------------------------------------
%           optk =     optimum number of clusters by K&L criterion.
%           centr =    [optk x p] matrix of cluster centroids.
%           clst =     [n x maxk] matrix of cluster memberships, in which column i 
%                        gives the cluster memberships for the case of i clusters.
%           Cg =       [maxk x 1] vector of Cg criterion values for k = 1:maxk
%                        clusters.
%

% Hardy,A. 1994. An examination of procedures for determining the number of
%   clusters in a data set.  In: Diday et al. (eds.), New Approaches in
%   Classification and Data Analysis, pp. 178-185.  Springer-Verlag.
%
% Krzanowski,WJ & YT Lai. 1988. A criterion for determining the number of
%   groups in a data set using sum-of-squares clustering.  Biometrics 44:23-34.

% RE Strauss, 8/17/98
%   9/7/99 - miscellaneous changes for Matlab v5.
%   10/22/02 - update documentation; add 'doplot' option.

function [optk,centr,clst,Cg] = kmeangrp(X,maxk,restarts,doplot)
  if (nargin < 2) maxk = []; end;
  if (nargin < 3) restarts = []; end;
  if (nargin < 4) doplot = []; end;

  if (isempty(restarts))
    restarts = 0;
  end;
  if (isempty(maxk))
    maxk = 6;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  [n,p] = size(X);

  cum_centr = zeros(sum(2:(maxk+1)),p);
  cum_clst =  zeros(n,maxk);
  c = 1;

  C = zeros(maxk+1,1);
  C(1) = sum(diag(cov(X)));

  for k = 2:(maxk+1)
    [centr,clst] = kmeans(X,k,restarts);
    cum_centr(c:(c+k-1),:) = centr;
    cum_clst(:,k-1) = clst;
    c = c+k;

    meancent = X - centr(clst,:);
    C(k) = k^(2/p) * sum(diag(cov(meancent)));
  end;

  diff = zeros(maxk+1,1);
  diff(2:maxk+1) = C(1:(maxk)) - C(2:maxk+1);
  Cg = abs(diff(2:maxk) ./ diff(3:maxk+1));

  [maxCg,optk] = max(Cg);
  optk = optk + 1;

  c = sum(1:(optk-1));
  centr = cum_centr(c:(c+optk-1),:);
  clst = [ones(n,1) cum_clst(:,1:maxk-1)];

  if (maxk > 2 & doplot)
    figure;
    plot([2:maxk]',Cg(1:maxk-1),'k');
    putbnd([2:maxk]',Cg(1:maxk-1));
    puttick([2:maxk]);
    putxlab('Number of clusters');
    putylab('Cg criterion');
  end;

  return;
