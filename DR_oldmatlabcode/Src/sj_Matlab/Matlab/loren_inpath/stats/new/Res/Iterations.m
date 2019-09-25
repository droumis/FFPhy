% ITERATIONS: Estimates the minimum number of iterations required to estimate a 
%             two-sided p-value to within delta units of the true one with  
%             (1-alpha) percent confidence, assuming that the iterated samples 
%             are independent (Agresti, Wackerly & Boyett 1979).
%             Input arguments can be scalars or matrices.
%
%     Usage: iter = iterations({delta},{ci_level})
%
%             delta =     confidence regions about p [default=0.01].
%             ci_level =  confidence levels [default=95].
%             -------------------------------------------------------------------
%             iter =      minimum number of iterations (as integer).
%

% Agresti, A, D Wackerly & JM Boyett. 1979. Exact conditional tests for cross-
%   classification: approximation of attained significance levels.
%   Psychometrika 44:75-83.

% RE Strauss, 5/19/00

function iter = iterations(delta,ci_level)
  if (nargin < 1) delta = []; end;
  if (nargin < 2) ci_level = []; end;

  if (isempty(delta))
    delta = 0.01;
  end;
  if (isempty(ci_level))
    ci_level = 0.95;
  end;
  if (any(ci_level>1))
    ci_level = ci_level./100;
  end;
  ci_level = 1-((1-ci_level)./2);       % Make two-sided

  ds = isscalar(delta);
  dc = isscalar(ci_level);

  if (ds & ~dc)
    delta = delta * ones(size(ci_level));
  elseif (~ds & dc)
    ci_level = ci_level * ones(size(delta));
  elseif (~ds & ~dc)
    if (size(delta) ~= size(ci_level))
      error('  ITERATIONS: input matrices not compatible.');
    end;
  end;

  z = norminv(ci_level);
  iter = ceil((z./(2*delta)).^2);

  return;
