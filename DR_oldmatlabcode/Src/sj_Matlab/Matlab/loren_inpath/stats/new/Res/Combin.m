% COMBIN: Calculates the number of possible combinations of N objects
%         taken R at a time.  If N and R are compatible vectors, a
%         corresponding vector is returned based on element-wise operations.
%         If N is a scalar and R a vector, N is expanded to the same size as R.
%
%     Usage: c = combin(N,R)
%

% RE Strauss, 1/13/96
%   11/4/98 - for scalar N and vector R, expand N to vector.

function c = combin(N,R)
  if (nargin<2)
    error('  COMBIN: Too few input arguments');
  end;

  if (min(size(N))~=1 | min(size(R))~=1)
    error('  COMBIN: vector input only');
  end;

  nsize = length(N);
  rsize = length(R);

  if (nsize==1 & rsize>1)     % Expand N if a scalar
    N = N*ones(size(R));
    nsize = length(N);
  end;

  if (nsize ~= rsize)
    error('  COMBIN: Argument vectors incompatible');
  end;

  c = N;
  for i=1:nsize
    if (N(i)>=R(i))
      c(i) = prod(1:N(i))/(prod(1:(N(i)-R(i)))*prod(1:R(i)));
    else
      c(i) = 0;
    end;
  end;
  return;
