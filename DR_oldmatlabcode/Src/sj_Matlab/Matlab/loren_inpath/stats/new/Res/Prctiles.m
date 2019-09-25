% PRCTILES: Finds the data values (possibly interpolated) corresponding to the
%           specified percentiles, by column.  Uses linear interpolation between
%           empirical values.  Allows for missing data.
%
%     Usage: prc = prctiles(X,prcvals)
%
%           X =       [n x p] data matrix.
%           prcvals = vector (length v) of percentile values to be estimated.
%           -----------------------------------------------------------------
%           prc =     [v x p] matrix of percentiles.
%

% RE Strauss, 1/11/99
%   6/14/99 - allow for missing data.

function prc = prctiles(X,prcvals)
  [r,c] = size(prcvals);
  if (min([r,c])>1)
    error('PRCTILES: prcvals must be a vector');
  end;
  v = max([r,c]);

  if (any(prcvals<0 | prcvals>100))
    error('PRCTILES: prcvals must be within range 0-100');
  end;

  [r,c] = size(X);
  if (r==1 & c>1)                   % If X is a row vector,
    X = X';                           % Transpose
    [r,c] = size(X);
  end;

  xx = sort(X);

  if (any(~finite(xx(:))))            % If any nonfinite values in matrix
    prc = zeros(v,c);
    for ic = 1:c                        % Cycle thru columns
      xxc = xx(:,ic); 
      i = find(~finite(xxc));
      if (~isempty(i))                  % Remove any missing values
        xxc(i) = [];
      end;
      rx = length(xxc);
      if (rx>0)
        q = [0 100*(0.5:rx-0.5)./rx 100];
        xxc = [min(xxc); xxc; max(xxc)];
        prc(:,ic) = interp1(q,xxc,prcvals);
      else
        prc(:,ic) = NaN*ones(v,1);
      end;
    end;
  else                              % Else handle whole matrix simultaneously
    q = [0 100*(0.5:r-0.5)./r 100];
    xx = [min(X); xx; max(X)];
    prc = interp1(q,xx,prcvals);
  end;

  return;


