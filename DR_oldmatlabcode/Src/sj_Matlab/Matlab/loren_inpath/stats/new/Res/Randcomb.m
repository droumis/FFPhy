% RANDCOMB: Produces a series of random combinations of N integers taken r at a 
%           time.
%
%     Usage: combs = randcomb(N,r,{ncombs},{nodups})
%
%           N =       number of objects from which sampled.
%           r =       number of objects sampled.
%           ncombs =  number of random combinations returned.  If not provided, 
%                       all are returned in a random sequence.  If 'ncombs' > 
%                       number of possible combinations, some are repeated.
%           nodups =  optional boolean flag indicating no random combinations 
%                       are to be duplicated (i.e., sampled without replacement)
%                       [default = false = 0].
%           --------------------------------------------------------------------
%           combs =   [ncombs x r] matrix of random combinations (rows).
%

% RE Strauss, 3/5/99

function combs = randcomb(N,r,ncomb,nodups)
  if (nargin < 3)
    ncomb = [];
  end;
  if (nargin < 4)
    nodups = [];
  end;

  if (~isintegr([N r]) | min([N r])<0 | r>N)
    error('RANDCOMB: N and r must be positive integers, r <= N');
  end;

  ncombpossible = combin(N,r);
  if (isempty(ncomb))
    ncomb = ncombpossible;
  end;

  if (isempty(nodups))
    nodups = 0;
  end;

  if (ncomb == ncombpossible)       % If want all possible combinations,
    combs = combvals(N,r);            % Get all possible, in order
    i = randperm(ncomb);              % Randomly permute
    combs = combs(i,:);

  elseif (ncomb > ncombpossible)    % If want more than all possible,
    combs = zeros(ncomb,r);
    c = combvals(N,r);                % Get all possible, in order
    i = randperm(ncombpossible);
    c = c(i,:);                       % Randomly permute

    b = 1;
    e = ncombpossible;
    reps = floor(ncomb/ncombpossible);
    for j = 1:reps                    % Append sets of all possible
      combs(b:e,:) = c;
      b = e+1;
      e = b+ncombpossible-1;
    end;

    mod = rem(ncomb,reps*ncombpossible);
    if (mod>0)
      combs(b:(b+mod-1),:) = c(i(1:mod),:);
    end;

    i = randperm(ncomb);              % Randomly permute entire set
    combs = combs(i,:);
  else                              % If want fewer than all possible,
    max_feasible = 5000;   

    if (ncombpossible < max_feasible) % If feasible to store all possible,
      combs = combvals(N,r);            % Get all possible, in order
      i = randperm(ncombpossible);      % Randomly permute
      i = i(1:ncomb);                   % Save subset
      combs = combs(i,:);
      
    else
      combs = zeros(ncomb,r);           % Allocate matrix
      if (nodups)                       % Prepare to check for duplicates
        base = N+1;
        decval = zeros(ncomb,1);
      end;

      for i = 1:ncomb                   % Generate random combinations
        p = randperm(N);
        p = sort(p(1:r));
        if (nodups)
          d = p * [base.^(r-1:-1:0)]';    % Convert to decimal identifier
          while (sum(d==decval)>0)
            p = randperm(N);
            p = sort(p(1:r));
            d = p * [base.^(r-1:-1:0)]';
          end;
          decval(i) = d;
        end;
        combs(i,:) = p;
      end;
    end;
  end;

  return;
