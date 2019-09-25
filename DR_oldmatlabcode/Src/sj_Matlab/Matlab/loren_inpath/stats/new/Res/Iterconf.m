% ITERCONF: Estimates confidence interval about a randomized p-value as a 
%           function of number of iterations, assuming that the iterated samples 
%           are independent (Agresti, Wackerly & Boyett 1979).
%           Input arguments can be scalars or matrices.
%
%     Usage: conf = iterconf(iter,{ci_level})
%
%             iter = number of randomization iterations.
%             ci_level =  confidence levels [default=95].
%             -------------------------------------------
%             conf = symmetric confidence bounds about p.


% Agresti, A, D Wackerly & JM Boyett. 1979. Exact conditional tests for cross-
%   classification: approximation of attained significance levels.
%   Psychometrika 44:75-83.

% RE Strauss, 5/19/00

function conf = iterconf(iter,ci_level)
  if (nargin < 2) ci_level = []; end;

  if (isempty(ci_level))
    ci_level = 0.95;
  end;
  if (any(ci_level>1))
    ci_level = ci_level./100;
  end;
  ci_level = 1-((1-ci_level)./2);       % Make two-sided

  di = isscalar(iter);
  dc = isscalar(ci_level);

  if (di & ~dc)
    iter = iter * ones(size(ci_level));
  elseif (~di & dc)
    ci_level = ci_level * ones(size(iter));
  elseif (~di & ~dc)
    if (size(iter) ~= size(ci_level))
      error('  ITERCONF: input matrices not compatible.');
    end;
  end;

  z = norminv(ci_level);
  conf = z./(2*sqrt(iter));

  return;


