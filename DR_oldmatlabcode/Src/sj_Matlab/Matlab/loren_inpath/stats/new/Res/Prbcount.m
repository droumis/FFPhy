% PRBCOUNT: Given a vector of probabilities and a total count, randomly or
%           exactly replaces the probabilities with counts summing to the total.
%           Useful for providing marginal totals, given the marginal probabilities.
%
%     Syntax: counts = prbcount(probs,N,{maxcount},{letzeros},{exact})
%
%           probs =    vector of cell probabilities.
%           N =        scalar indicating total count.
%           maxcount = maximum cell count [default = N].
%           letzeros = boolean flag indicating whether (=1) or not (=0) zero cell
%                       counts are to be permitted [default = 0].  If any
%                       cell probability is zero, letzero is set to zero.
%           exact =    boolean flag indicating that counts are to be as exact
%                       as possible (=1) rather than randomized (=0) [default = 0].
%           -----------------------------------------------------------------------
%           counts =  corresponding vector of cell counts.
%

% RE Strauss, 4/17/97
%   9/19/99 -  updated handling of default input arguments.
%   10/16/01 - added error message for small N.
%   12/10/01 - changed condition test for small N error message.
%   8/9/03 -   allow for zero probabilities.

function counts = prbcount(probs,N,maxcount,letzeros,exact)
  if (nargin < 3) maxcount = []; end;
  if (nargin < 4) letzeros = []; end;
  if (nargin < 5) exact = []; end;

  [r,c] = size(probs);
  if (min([r,c])>1)
    error('  PRBCOUNT: Probabilities must be in vector form');
  end;
  if (abs(sum(probs)-1) > 1e-6)
    error('  PRBCOUNT: Probabilities must sum to unity');
  end;
  if (max(size(N))>1)
    error('  PRBCOUNT: N must be a scalar');
  end;

  if (isempty(maxcount)) maxcount = N; end;
  if (isempty(letzeros)) letzeros = 0; end;
  if (isempty(exact))    exact = 0;    end;

  if (any(probs < 1e-6))
    letzeros = 1;
  end;

if (~letzeros & N<length(probs))
    error('   PRBCOUNT: N too small');
  end;

  cumprobs = cumsum(probs);

  if (exact)                            % Convert probabilities to exact counts
    counts = round(probs*N);
    while (sum(counts)~=N)
      dev = probs*N - counts;
      if (sum(counts) > N)
        [maxdev,i] = max(-dev);
        i = i(1);
        counts(i) = counts(i)-1;
      else
        [maxdev,i] = max(dev);
        i = i(1);
        counts(i) = counts(i)+1;
      end;
    end;

  else
    counts = zeros(size(probs));        % Convert probabilities to randomized counts
    randval = rand(N,1);
    for i = 1:N
      f = max(find(randval(i)>cumprobs))+1;
      if (isempty(f))
        f = 1;
      end;
      counts(f) = counts(f)+1;
    end;
  end;

  if (~letzeros)                        % Check for zero counts
    lenprobs = length(probs);
    while(any(counts<1))
      z = min(find(counts<1));            % If more than 1, choose the 1st
      f = z;
      while (f==z)                        % Pick a random cell based on probs
        f = max(find((rand*(lenprobs+1)/lenprobs)>cumprobs));
        if (isempty(f))
          f = 1;
        end;
      end;
      if (counts(f)>1)
        counts(f) = counts(f)-1;            % Move an observation
        counts(z) = counts(z)+1;
      end;
    end;
  end;

  while (any(counts>maxcount))        % Check for counts above max
    z = min(find(counts>maxcount));     % If more than 1, choose the first
    f = z;
    while (f==z)                        % Pick a random cell based on probs
      f = max(find(rand>cumprobs));
      if (isempty(f))
        f = 1;
      end;
    end;
    if (counts(f)<maxcount)
      counts(f) = counts(f)+1;            % Move an observation
      counts(z) = counts(z)-1;
    end;
  end;

  return;
