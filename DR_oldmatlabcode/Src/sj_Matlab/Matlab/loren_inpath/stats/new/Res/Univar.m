% UNIVAR: For each column in a data matrix, returns a vector
%         containing the mean, median, average deviation from the median,
%         standard deviation, variance, skewness, and kurtosis.
%         Ignores missing values (unlike the corresponding Matlab functions).
%         The unibased estimate of the standard deviation is from Winkler &
%         Hays, 1970).
%
%     Usage: stats = univar(X,{get_stats})
%
%         X =         [n x p] data matrix, analyzed by column.
%         get_stats = vector of boolean flags indicating which statistics are to be
%                       estimated [default = all].
%         ------------------------------------------------------------------
%         stats =     [s x p] matrix containing a column vector of
%                        statistics for each requested variable:
%                           1) mean
%                           2) median
%                           3) mean absolute deviation from the median
%                           4) standard deviation
%                           5) variance
%                           6) skewness
%                           7) kurtosis
%                           8) sample size (excluding non-finite data)
%

% RE Strauss
%   1/16/99 - added boolean vector of options
%   6/24/99 - allow for all missing data

function stats = univar(X,get_stats)
  if (~nargin) help univar; return; end;

  if (nargin < 2) get_stats = []; end;

  nstats = 8;
  if (isempty(get_stats))
    get_stats = ones(1,nstats);
  else
    if (size(get_stats,1)>1)          % Pad get_stats vector to full length
      get_stats = get_stats';
    end;
    if (length(get_stats) < nstats)
      get_stats = [get_stats, zeros(1,nstats-length(get_stats))];
    end;
  end;

  [N,P] = size(X);                    % Numbers of observations & variables
  if (N == 1)
    if (P > 1)                        % Transpose input row matrix
      X = X';
      [N,P] = size(X);
    end;
  end;

  for c = 1:P                         % Cycle thru columns
    x = X(:,c);
    b = ~finite(x);
    x(b) = [];                          % Remove missing data

    nx = length(x);
    sampsize(c) = nx;

    if (get_stats(1))                   % Statistics
      if (nx >= 1)
        meanx(c) = mean(x);
      else
        meanx(c) = NaN;
      end;
    end;
    if (get_stats(4) | get_stats(5))
      if (nx >= 2)
        varx(c) = var(x);
      else
        varx(c) = NaN;
      end;
    end;
    if (get_stats(4))
      if (nx >= 2)
        stdev(c) = sqrt(varx(c)) .* (1 + 1./((4.*sampsize(c))-4));
      else
        stdev(c) = NaN;
      end;
    end;
    if (get_stats(6))
      if (nx >= 3)
        skew(c) = skewness(x);
      else
        skew(c) = NaN;
      end;
    end;
    if (get_stats(7))
      if (nx >= 4)
        kurt(c) = kurtosis(x);
      else
        kurt(c) = NaN;
      end;
    end;
    if (get_stats(2) | get_stats(3))
      if (nx >= 1)
        med(c) = median(x);
      else
        med(c) = NaN;
      end;
    end;
    if (get_stats(3))
      if (nx >= 1)
        absdevs(c) = sum(abs(x-ones(sampsize(c),1)*med(c)))./N;
      else
        absdevs(c) = NaN;
      end;
    end;
  end;

  stats = [];                       % Prepare output matrix
  if (get_stats(1))
    stats = [stats; meanx];
  end;
  if (get_stats(2))
    stats = [stats; med];
  end;
  if (get_stats(3))
    stats = [stats; absdevs];
  end;
  if (get_stats(4))
    stats = [stats; stdev];
  end;
  if (get_stats(5))
    stats = [stats; varx];
  end;
  if (get_stats(6))
    stats = [stats; skew];
  end;
  if (get_stats(7))
    stats = [stats; kurt];
  end;
  if (get_stats(8))
    stats = [stats; sampsize];
  end;

  return;
