% KWSTAT: Calculates the Kruskal-Wallace statistic value, given either the data 
%         or the ranks of the data.
%
%     Usage: H = kwstat(x,g,are_ranks)
%
%         x =         vector (length n) of data or ranks of data.
%         g =         corresponding grouping vector.
%         are_ranks = boolean flag indicating that 'x' consists of ranks 
%                       [default = 0 = false].
%         --------------------------------------------------------------
%         H =         K-W statistic value.
%

% Sokal & Rohlf, 1981, pp. 430-432.

% RE Strauss, 10/27/99

function H = kwstat(x,g,are_ranks)
  if (nargin < 3) are_ranks = []; end;

  if (isempty(are_ranks))
    are_ranks = 0;
  end;

  g = g(:);
  if (are_ranks)
    r = x(:);
  else
    r = ranks(x(:));
  end;

  [ur,fr] = uniquef(r);                 % Unique ranks and their freqs
  [ug,n] = uniquef(g);                  % Groups and their sample sizes
  ngrps = length(ug);                   % Number of groups
  N = length(g);                        % Total sample size

  sumr = zeros(ngrps,1);                % Sums of ranks, by group
  for ig = 1:ngrps
    i = find(g==ug(ig));
    sumr(ig) = sum(r(i));
  end;

  sr = sum((sumr.^2)./n);               % Test-statistic value
  H = 12*sr./(N*(N+1)) - 3*(N+1);

  if (any(fr>1))                        % If there are ties,
    t = fr(fr>1);                         % Number of ties in each group of ties
    T = (t.^3) - t;
    D = 1 - sum(T)/(N.^3 - N);
    H = H/D;
  end;

  return;
