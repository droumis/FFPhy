% MODE: Heuristic for determining the mode of a distribution.
%       If data consist of integers: determines frequencies; if max frequencies 
%         are contiguous, returns the mean of the values having max frequencies, 
%         otherwise returns NaN.
%       If data are non-integer: bins data into numbers of bins varying from 
%         3:(N-3); for each binning yielding a single max frequency, saves the 
%         number of bins and the class midpoint; estimated mode is the mean 
%         class midpoint, weighted by the deviation of the number of bars from 
%         the 'best' number of bars based on Sturges' criterion.  If no binning 
%         yields a single max frequency, returns NaN.
%
%     Usage: m = mode(X)
%
%         X = [n x p] data matrix
%         ---------------------------------------------
%         m = [1 x p] vector of modes for columns of X.
%

% RE Strauss, 2/25/01 (replaced previous version of 5/8/95)

function m = mode(X)
  if (isvector(X))
    X = X(:);
  end;
  [N,P] = size(X);

  m = NaN*ones(1,P);

  for col = 1:P
    x = X(:,col);
    if (isintegr(x))                  % Integer data
      [u,f] = uniquef(x);
      i = find(f==max(f));
      leni = length(i);
      if (i==1)
        m(col) = u(i);
      else
        d = i(2:leni)-i(1:(leni-1));
        if (all(d==ones(1,leni-1)))
          m(col) = mean(u(i));
        end;
      end;

    else                              % Continuous data
      optbars = ceil((16.0*N)^(1/3)+0.5);
      bars = [];
      modeguess = [];

      for nbins = 3:(N-3)
        [f,midpt] = hist(x,nbins);
        i = find(f==max(f));
        if (length(i)==1)
          bars = [bars; nbins];
          modeguess = [modeguess; midpt(i)];
        end;
      end;
      if (~isempty(modeguess))
        m(col) = meanwt(modeguess,optbars-abs(bars-optbars));
      end;
    end;
  end;
