% CENTDIST: Given a vector x, returns the values of x delimiting the central 
%           alpha-percent of the distribution (interpolated if necessary).
%
%     Usage: lims = centdist(X,{alpha})
%
%           X =     [n x p] matrix in which columns represent distributions.
%           alpha = percentage of distribution to be captured [default = 95].
%           -----------------------------------------------------------------
%           lims =  [2 x p] matrix specifying the lower (row 1) and upper (row 2) 
%                     values capturing the central portion of the distributions.
%

function lims = centdist(X,alpha)
  if (nargin < 2)
    alpha = [];
  end;
  if (isempty(alpha))
    alpha = 0.95;
  else
    if (alpha > 1)
      alpha = 0.01*alpha;
    end;
  end;


  [n,p] = size(X);
  if (n == 1)
    X = X';
    [n,p] = size(X);
  end;
  lims = zeros(2,p);                  % Allocate return vector

  X = sort(X);



  for v = 1:p;                        % Cycle thru variables
    D = X(:,v);                         % Extract distribution

    if (var(D)>0)                         % If distribution is not invariant,
      L = (1-alpha)/2;                    % Percentiles
      U = 1-L;

      indx = n * L;                       % Find lower limit
      low =  min(max(floor(indx),1),n);
      high = max(min(ceil(indx),n),1);

      if (high-low > 0)
        delta = (indx - low)/(high - low);
        lim(1,v) = D(low) + delta*(D(high)-D(low));
      else
        lim(1,v) = D(low);
      end;

      indx = n * U;                       % Find upper limit
      low =  min(max(floor(indx),1),n);
      high = max(min(ceil(indx),n),1);

      if (high-low > 0)
        delta = (indx - low)/(high - low);
        lims(2,v) = D(low) + delta*(D(high)-D(low));
      else
        lims(2,v) = D(low);
      end;
    else
        lims(1,v) = X(1,v);
        lims(2,v) = X(1,v);
    end;
  end;

  return;



