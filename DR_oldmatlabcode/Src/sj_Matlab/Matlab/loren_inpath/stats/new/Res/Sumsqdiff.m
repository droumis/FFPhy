% SUMSQDIFF:  Calculates distances among rows (or cols) as the sum of squared 
%             differences among cols (or rows), optionally bootstrapped.
%
%     Usage: [dist,CI_dist] = sumsqdiff(X,{among_cols},{iter},{CI_level})
%
%             X =           [r x c] matrix.
%             among_cols =  optional boolean flag indicating, if true, that 
%                             distances are to be calculated among columns 
%                             instead of among rows [default=0].
%             iter =        optional number of bootstrap iterations [default=0].
%             CI_level  =   optional level of bootstrapped confidence intervals 
%                             [default = 95]
%             ------------------------------------------------------------------
%             dist =        [r x r] matrix of pairwise distances among rows 
%                             (or [c x c] matrix among cols).
%             CI_dist =     [r x r] matrix of bootstrapped confidence intervals, 
%                             with lower limits below the diagonal and upper 
%                             limits above the diagonal.
%

% RE Strauss, 5/16/00

function [dist,CI_dist] = sumsqdiff(X,among_cols,iter,CI_level)
  if (nargin < 2) among_cols = []; end;
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  if (isempty(among_cols))
    among_cols = 0;
  end;
  if (isempty(iter) | nargout < 2)
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  end;
  if (CI_level > 1)
    CI_level = CI_level/100;
  end;

  if (among_cols)
    X = X';
  end;

  [r,c] = size(X);
  dist = zeros(r,r);
  CI_dist = [];

  for i = 1:(r-1)                         % Distance matrix
    for j = (i+1):r
      d = X(i,:) - X(j,:);
      d = d*d';
      dist(i,j) = d;
      dist(j,i) = d;
    end;
  end;

  if (iter)                               % Bootstrapped confidence intervals
    ci = bootstrp('sumsqdifff',[1 0 0 1],iter,1-CI_level,X');
    CI_dist = trisqmat(ci');
  end;

  return;
