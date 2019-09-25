% DATOGRPS: Given a data matrix and vector of group identifiers, performs a discriminant 
%           or size-free discriminant analysis and passes the mean scores (by group) to the 
%           kmeansgrp function to estimate the optimal number of groups.
%
%     Usage: [optk,centr,clst,Cg] = 
%                 datogrps(X,grps,{kindanal},{nfact},{maxk},{restarts},{kindregr})
%
%           X -         [n x p] data matrix (obs x vars).
%           grps -      row or column vector of group identifiers.
%           kindanal -  kind of multivariate analysis:
%                         0 = PCA (cov matrix)
%                         1 = DFA [default]
%                         2 = size-free DFA
%           nfact -     optional number of leading discriminant functions for
%                         which scores are desired [default=3].
%           maxk =      maximum number of clusters fitted [default = number of groups-1].
%           restarts =  number of restarts, after finding a new minimum sse,
%                         needed to end search for each value of k [default = 0].
%           kindregr =  regression model used (0=major axis [default], 1=predictive) for 
%                         the size-free discriminant analysis procedure, if used.
%           ----------------------------------------------------------------------------
%           optk =      optimum number of clusters by K&L criterion.
%           centr =     [optk x p] matrix of cluster centroids.
%           clst =      [n x 1] vector of cluster memberships.
%           Cg =        [maxk x 1] vector of Cg criterion values for k = 1:maxk
%                         clusters.
%

% RE Strauss, 7/21/98
%   11/29/99 - changed calling sequence.

function [optk,centr,clst,Cg] = datogrps(X,grps,kindanal,nfact,maxk,restarts,kindregr)

  if (nargin < 3) kindanal = []; end;
  if (nargin < 4) nfact = []; end;
  if (nargin < 5) maxk = []; end;
  if (nargin < 6) restarts = []; end;
  if (nargin < 7) kindregr = []; end;

  g = uniquef(grps,1);                  % Determine groups
  ngrps = length(g);

  if (isempty(kindanal))
    kindanal = 1;
  end;
  if (isempty(nfact))
    nfact = 3;
  end;
  if (isempty(maxk))
    maxk = ngrps-1;
  end;
  
  if (kindanal == 0)                    % Ordination
    [load,percvar,scores] = pcacov(X,nfact);
  elseif (kindanal == 1)
    [load,percvar,scores] = discrim(X,grps,nfact);
  elseif (kindanal == 2)
    [load,percvar,scores] = sizefree(X,grps,nfact,[],[],[],kindregr);
  else
    error('  DATOGRPS: invalid kind of analysis');
  end;

  scores = zscore(scores);              % Standardize scores by factor
  for f = 1:nfact                       % Scale factors by percent variance
    scores(:,f) = scores(:,f) * percvar(f);
  end;

  centr = zeros(ngrps,nfact);           % Centroids by group
  for i = 1:ngrps
    centr(i,:) = mean(scores(grps==g(i),:));
  end;
  
  [optk,centr,clst,Cg] = kmeangrp(centr,maxk,restarts);   % Find optimal number of groups
  grp_cluster = [g clst]

  return;

