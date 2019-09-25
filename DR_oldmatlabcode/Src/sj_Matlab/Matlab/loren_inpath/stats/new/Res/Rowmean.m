% ROWMEAN: Returns column vector of row means of a matrix, optionally ignoring 
%          missing data.  Returns a scalar for an row vector, and the original 
%          input for a column vector.
%
%     Usage: r = rowmean(X,{miss})
%
%         X =     [n x p] matrix.
%         miss =  optional boolean flag indicating that missing (non-finite) 
%                   values are to be ignored [default = 0].  If true, returns 
%                   NaN values only for rows having all non-finite values.
%         -------------------------------------------------------------------
%         r =     column vector of row sums.
%

% RE Strauss, 6/4/01 (modified from rowsum.m)

function r = rowmean(X,miss)
  if (~nargin) help rowmean; return; end;

  if (nargin < 2) miss = []; end;

  if (isempty(miss))
    miss = 0;
  end;

  if (miss)
    [isvectX,ncells,iscol] = isvector(X);
    if (isvectX)
      if (iscol)
        r = X;
      else
        X(~isfinite(X)) = [];
        r = mean(X);
      end;
    else
      [n,p] = size(X);
      r = zeros(n,1);
      for ir = 1:n
        x = X(ir,:);
        x(~isfinite(x)) = [];
        if (isempty(x))
          r(ir) = NaN;
        else
          r(ir) = mean(x);
        end;
      end;
    end;
  else
    r = mean(X')';
  end;

  return;
