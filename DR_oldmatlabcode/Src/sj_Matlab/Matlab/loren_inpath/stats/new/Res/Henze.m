% HENZE: Applies the Henze (1988) distribution-free test for the 
%        inequality of two or more multivariate distributions, based on 
%        nearest-neighbor coincidences.
%        Randomized significance level is evaluated by permuting group labels.
%
%     Syntax: [prob,coincidences] = henze(X,grps,iter)
%
%           X -    [N x p] data matrix.
%           grps - [N x 1] vector of group labels.
%           iter - number of permutation iterations [default = 5000].
%           -----------------------------------------------------------
%           prob -         significance level.
%           coincidences - observed number of multivariate runs.
%

% Good, PG. 1994. Permutation Tests, p. 72. Springer Series in Statistics,
%   Springer-Verlag.

% RE Strauss, 7/4/96
%   11/29/99 - changed calling sequence.

function [prob,coincidences] = henze(X,grps,iter)
  if (nargin < 3) iter = []; end;

  if (isempty(iter))
    iter = 5000;
  end;

  N = length(grps);

  dist = eucl(X);                     % Euclidean dists among pts
  dist = putdiag(dist,1e12);          % Put high values on diagonal

  [m,NN] = min(dist);                 % Find nearest neighbor of each pt
  coincidences = sum(grps==grps(NN)); % Number of times NN is in same grp

  prob = 0;
  incr = 1/iter;

  for it = 1:iter                     % Randomized significance level
    grp = grps(randperm(N));            % Randomly permute labels
    c = sum(grp==grp(NN));              % Number of times NN is in same grp

    if (c >= coincidences)              % Increment tail probability
      prob = prob + incr;
    end;
  end;

  return;

